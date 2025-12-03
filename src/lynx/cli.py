import argparse
import sys
import os
import gzip
import re
import subprocess
import logging

from collections import defaultdict
from itertools import chain
from statistics import median
from functools import partial, partialmethod

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from tqdm.contrib.logging import logging_redirect_tqdm

from rich_argparse import ArgumentDefaultsRichHelpFormatter
from . import __version__

logging.basicConfig(
    level='INFO',
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

if sys.stderr.isatty():
    logging.addLevelName(logging.INFO, f'\033[1m{logging.getLevelName(logging.INFO)}\033[1;0m')
    logging.addLevelName(logging.WARNING, f'\033[1m\x1b[33;20m{logging.getLevelName(logging.WARNING)}\033[1;0m')
    logging.addLevelName(logging.CRITICAL, f'\033[1m\x1b[31;20m{logging.getLevelName(logging.CRITICAL)}\033[1;0m')
else:
    tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

log = logging.getLogger(__name__)

ArgumentDefaultsRichHelpFormatter.styles['argparse.prog'] = 'default'
ArgumentDefaultsRichHelpFormatter.styles['argparse.default'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.metavar'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.groups'] = '#C19A6B'
ArgumentDefaultsRichHelpFormatter.styles['argparse.args'] = 'default'


def cli():
    parser = argparse.ArgumentParser(
        prog='lynx',
        description=f'Lynx v{__version__}: lightweight profiling of antibiotic resistance genes from short-read metagenomes',
        formatter_class=ArgumentDefaultsRichHelpFormatter
    )

    parser.add_argument(
        dest='files',
        nargs='+',
        metavar='file',
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzip optional <*.gz>.')

    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-d',
        '--db',
        metavar='FILE',
        required=True,
        help='Database file.')

    required.add_argument(
        '-o',
        '--out',
        metavar='DIR',
        required=True,
        help='Output folder.')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        '--force',
        action='store_true',
        help='Force overwriting.')

    optional.add_argument(
        '--single',
        action='store_true',
        help='Files are single-end. If not given then merge forward|reward files with <_(1|2)>, <_(R1|R2)> or <_(fwd|rev)>.'
    )

    optional.add_argument(
        '-i',
        '--id',
        metavar='FLOAT',
        type=float,
        default=90,
        help='Min. identity in percentage.')

    optional.add_argument(
        '-c',
        '--cov',
        metavar='FLOAT',
        type=float,
        default=70,
        help='Min. query cover in percentage.')

    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help='Number of threads.')

    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.set_defaults(func=run)

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.func(**{key: val for key, val in vars(args).items() if key != 'func'})


def parse(file, id, cov):
    a = {'l15e', 's19e', 'l4e', 's24e', 'l10e', 'l44e', 's3ae', 'l11'}
    b = {'s7', 'l2', 'l11', 'l13', 's2', 's16', 'l27', 'l20'}
    adict = {x: 0 for x in a}
    bdict = {x: 0 for x in b}
    xdict = defaultdict(float)
    ydict = defaultdict(float)
    zdict = defaultdict(float)

    cov /= 100
    with gzip.open(file, 'rt') as f:
        cur_q, block = None, []
        for line in chain(f, [None]):
            ls = None if line is None else line.split()
            qseqid = None if ls is None else ls[0]
            if cur_q is None and qseqid is not None:
                cur_q = qseqid

            if qseqid != cur_q:
                if block:
                    max_length = 0
                    scov_at_max = []
                    for r in block:
                        sseqid = r[1]; pident = float(r[2]); length = int(r[3]); qlen = int(r[4]); scov = float(r[5])
                        if sseqid[:4] == 'SARG':
                            if pident < id or length < max(25, min(50, (qlen / 3)) * cov):
                                continue
                        else:
                            if pident < 70 or length < max(25, min(50, (qlen / 3)) * 0.7):
                                continue

                        if length > max_length:
                            sseqid_at_max = sseqid
                            scov_at_max = [scov]
                            max_length = length

                        elif length == max_length:
                            scov_at_max.append(scov)

                    if scov_at_max:
                        s = median(scov_at_max) / 100
                        if 'SARG' in sseqid_at_max:
                            xdict[sseqid_at_max.split('|', 1)[-1]] += s
                            ydict[sseqid_at_max.split('|', 1)[-1].rsplit('|', 1)[0]] += s
                            zdict[sseqid_at_max.split('|', 1)[-1].rsplit('|', 2)[0]] += s
                        else:
                            m = sseqid_at_max.split('-')[0]
                            if 'bacteria' in sseqid_at_max and m in b:
                                bdict[m] += s
                            if 'archaea' in sseqid_at_max and m in a:
                                adict[m] += s

                block, cur_q = [], qseqid
                if qseqid is None:
                    break
            if ls:
                block.append(ls)

    genome = (sum(bdict.values()) / 8) + (sum(adict.values()) / 8)
    return os.path.basename(file).split('.tab.gz')[0], genome, xdict, ydict, zdict


def check(files, db, out, single, force):
    os.makedirs(out, exist_ok=True)
    for file in files + [db]:
        if not os.path.isfile(file):
            log.critical(f'File <{file}> does not exist.')
            sys.exit(2)

    samples = [re.sub('.gz$', '', os.path.basename(file)) for file in files]
    extensions = {sample.split('.')[-1] for sample in samples}
    if len(extensions) != 1 or not {'fasta', 'fa', 'fastq', 'fq'} & extensions:
        log.critical('Input files need to end with <fa|fq|fasta|fastq>.')
        sys.exit(2)
    extension = extensions.pop()

    if not single:
        samples = [re.sub(rf'(_(1|2|R1|R2|fwd|rev))?.{extension}$', '', sample) for sample in samples]
    else:
        samples = [re.sub(rf'.{extension}$', '', sample) for sample in samples]

    items = defaultdict(list)
    for sample, file in sorted(zip(samples, files)):
        items[sample].append(file)

    for sample, files in items.items():
        if len(files) > 2:
            log.critical(f'Sample <{sample}> has more than two input files. Check file format.')
            sys.exit(2)

    if len({len(file) for file in items.values()}) != 1:
        log.warning('Files are mixed with single-/paired-end. Check whether <--single> is needed.')

    return items


def run(files, db, out, single, force, id, cov, threads):
    with logging_redirect_tqdm():
        items = check(files, db, out, single, force)
        for sample, file in tqdm(items.items(), desc='Processing', leave=False):
            if os.path.isfile(f'{out}/{sample}.tab.gz') and not force:
                log.info(f'File <{out}/{sample}.tab.gz> exists, skip. Use <--force> for overwriting.')
            else:
                with open(f'{out}/{sample}.tab.gz', 'wb') as f:
                    for query in file:
                        p1 = subprocess.Popen([
                            'diamond', 'blastx',
                            '--query', query,
                            '--db', db,
                            '--threads', str(threads),
                            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'scovhsp', 'evalue', 'bitscore',
                            '--top', '0',
                            '--id', str(min(70, id)),
                        ], stdout=subprocess.PIPE)

                        p2 = subprocess.Popen(['pigz', '-c', '--best', '-p', str(threads)], stdin=p1.stdout, stdout=f)
                        p1.stdout.close()
                        p2.wait()
                        p1.wait()

    worker = partial(parse, id=id, cov=cov)
    with open(f'{out}/sequence.tsv', 'w') as f1, open(f'{out}/subtype.tsv', 'w') as f2, open(f'{out}/type.tsv', 'w') as f3:
        f1.write('sample\tsequence\tcopy\tgenome\tabundance\n')
        f2.write('sample\tsubtype\tcopy\tgenome\tabundance\n')
        f3.write('sample\ttype\tcopy\tgenome\tabundance\n')

        for chunk in process_map(worker, [f'{out}/{file}.tab.gz' for file in items.keys()], max_workers=threads, chunksize=1, desc='Parsing', leave=False):
            sample, genome, xdict, ydict, zdict = chunk
            for k, v in sorted(xdict.items()):
                f1.write(f"{sample}\t{k}\t{v}\t{genome}\t{(v / genome) if genome else 'n/a'}\n")
            for k, v in sorted(ydict.items()):
                f2.write(f"{sample}\t{k}\t{v}\t{genome}\t{(v / genome) if genome else 'n/a'}\n")
            for k, v in sorted(zdict.items()):
                f3.write(f"{sample}\t{k}\t{v}\t{genome}\t{(v / genome) if genome else 'n/a'}\n")

if __name__ == '__main__':
    cli()

import argparse
import sys
import os
import gzip
import re
import subprocess
import logging

from collections import defaultdict, Counter
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
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzipped optional <*.gz>.')

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
        help='Files are single-end. If not given then merge forward|reverse files with <_(1|2)>, <_(R1|R2)> or <_(fwd|rev)>.')

    optional.add_argument(
        '-l',
        '--level',
        metavar='STR',
        type=str,
        default='subtype',
        help='Aggregate ARG copies and abundance by <type|subtype|sequence>.')

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


def parse(file, id, cov, a, b, level=2):
    try:
        acopy = {x: 0 for x in a}
        bcopy = {x: 0 for x in b}
        scopy = defaultdict(float)

        cov /= 100
        with gzip.open(file, 'rt') as f:
            best = None
            old_qseqid = None
            for line in chain(f, [None]):
                ls = None if line is None else line.split()
                qseqid = None if ls is None else ls[0]

                if old_qseqid is None and qseqid is not None:
                    old_qseqid = qseqid

                if qseqid != old_qseqid:
                    if best is not None:
                        sseqid_at_max, max_length, scov_at_max = best
                        if scov_at_max:
                            s = median(scov_at_max) / 100
                            if sseqid_at_max.startswith('SARG'):
                                group = sseqid_at_max.split('|', 3)[1: level + 1]
                                group[0] = group[0].split('@', 1)[0]
                                scopy[tuple(group)] += s
                            else:
                                m = sseqid_at_max.split('-')[0]
                                if 'bacteria' in sseqid_at_max and m in b:
                                    bcopy[m] += s
                                if 'archaea' in sseqid_at_max and m in a:
                                    acopy[m] += s

                    if qseqid is None:
                        break

                    best = None
                    old_qseqid = qseqid

                if ls:
                    sseqid = ls[1]; pident = float(ls[2]); length = int(ls[3]); qlen = int(ls[4]); scov = float(ls[5])
                    if sseqid.startswith('SARG'):
                        if pident < id or length < max(25, min(50, qlen / 3) * cov):
                            continue
                    else:
                        if pident < 70 or length < max(25, min(50, qlen / 3) * 0.7):
                            continue

                    if best is None or length > best[1]:
                        best = (sseqid, length, [scov])
                    elif length == best[1]:
                        best[2].append(scov)

        genome = (sum(bcopy.values()) / 8) + (sum(acopy.values()) / 8)
        return os.path.basename(file).split('.out.gz')[0], genome, scopy

    except Exception as e:
        return file, None


def check(files, db, out, single, force, level):
    levels = ('type', 'subtype', 'sequence')
    if level not in levels:
        log.critical(f"Invalid --level <{level}>. Choose from <{'|'.join(levels)}>.")
        sys.exit(2)

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

    return items, levels[:levels.index(level) + 1]


def run(files, db, out, single, force, id, cov, threads, level='subtype'):
    threads = max(1, threads or 1)
    log.info('Aligning reads ...')
    with logging_redirect_tqdm():
        items, levels = check(files, db, out, single, force, level)
        for sample, file in tqdm(items.items(), leave=False):
            cache = f'{out}/{sample}.out.gz'
            if os.path.isfile(cache) and not force:
                log.warning(f'File <{cache}> exists, skipped. Use <--force> for overwriting.')
            else:
                tmp = f'{cache}.tmp'
                if os.path.exists(tmp):
                    os.remove(tmp)
                with open(tmp, 'xb') as f:
                    for query in file:
                        with subprocess.Popen([
                            'diamond', 'blastx',
                            '--query', query,
                            '--db', db,
                            '--threads', str(threads),
                            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'scovhsp', 'evalue', 'bitscore',
                            '--top', '0',
                            '--id', str(min(70, id)),
                            '--quiet',
                        ], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as p1:
                            p2 = subprocess.run(['pigz', '-c', '--best', '-p', str(threads)], stdin=p1.stdout, stdout=f)

                        if p1.returncode or p2.returncode:
                            os.remove(f.name)
                            log.critical(f'Failed to process <{query}>.')
                            sys.exit(2)

                os.replace(tmp, cache)

    log.info('Parsing outputs ...')
    acnt = Counter()
    bcnt = Counter()

    with open(db) as f:
        for line in f:
            if '>' in line and 'SARG' not in line:
                ls = line[1:-1].split('-')
                if ls[-2] == 'archaea':
                    acnt[ls[0]] += 1
                elif ls[-2] == 'bacteria':
                    bcnt[ls[0]] += 1

    a = {x[0] for x in acnt.most_common(8)}
    b = {x[0] for x in bcnt.most_common(8)}

    worker = partial(parse, id=id, cov=cov, a=a, b=b, level=len(levels))
    with open(f'{out}/output.abundance.tsv', 'w') as f1, open(f'{out}/output.genome.tsv', 'w') as f2:
        f1.write('\t'.join(('sample', *levels, 'copy', 'abundance')) + '\n')
        f2.write('sample\tgenome\n')

        for chunk in process_map(worker, [f'{out}/{file}.out.gz' for file in items.keys()], max_workers=threads, chunksize=1, leave=False):
            if chunk[-1] is None:
                log.warning(f"File <{chunk[0]}> failed to parse, skipped. Consider <--force> for overwriting.")
            else:
                sample, genome, scopy = chunk
                f2.write(f'{sample}\t{genome}\n')
                for group, copies in sorted(scopy.items()):
                    abundance = copies / genome if genome else 'n/a'
                    f1.write('\t'.join((sample, *group, str(copies), str(abundance))) + '\n')

    log.info('Finished.')

if __name__ == '__main__':
    cli()

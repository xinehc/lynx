# Lynx
Lynx: lightweight profiling of antibiotic resistance genes from short-read metagenomes

## Quick start
### Installation

<a href="https://anaconda.org/xinehc/lynx"> <img src="https://anaconda.org/xinehc/lynx/badges/version.svg" /> </a>
<a href="https://anaconda.org/xinehc/lynx"> <img src="https://anaconda.org/xinehc/lynx/badges/latest_release_date.svg" /> </a>

```bash
conda install -c bioconda -c conda-forge diamond pigz xinehc::lynx
```

### Database setup
> [!NOTE]
> In early versions of the ARG database, ARGs were uniquely identifiable only by type + subtype, not by subtype alone, because subtypes could be duplicated, e.g., *streptothricin satA* and *multidrug satA*. This has been fixed in databases released on or after 2026-05-20 by renaming certain subtypes to avoid naming collisions.

Download the latest database from [Zenodo](https://doi.org/10.5281/zenodo.17536666):
```bash
wget -qN --show-progress https://zenodo.org/records/20314770/files/db.fa
```

### Run Lynx
Download an example file containing 1,000,000 paired-end short reads:
```bash
wget -qN --show-progress https://zenodo.org/records/15681353/files/example.tar.gz
tar -xvf example.tar.gz
```

> [!TIP]
> If you have single-end reads, use `--single`; otherwise files with suffixes `<_(1|2)>`, `<_(R1|R2)>`, or `<_(fwd|rev)>` will be merged.

Run Lynx:
```bash
lynx -d db.fa -o . example/*.fasta
```

The `output.genome.tsv` file contains estimated genome copies obtained by mapping reads to eight bacterial/archaeal marker genes:

```text
sample	genome
sample	61.96256249999998
```

The `output.abundance.tsv` file contains estimated ARG copies and the associated abundances in copies per genome (i.e., ARG copies divided by genome copies):
```text
sample	type		subtype	 copy				 abundance
...
sample	sulfonamide	sul1	 4.16				 0.06713731376103113
sample	sulfonamide	sul2	 169.33200000000116	 2.7328114456209134
sample	sulfonamide	sul3	 26.763499999999997	 0.4319301675104222
...
```

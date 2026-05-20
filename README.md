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
<a href="https://doi.org/10.5281/zenodo.17536666"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17536666.svg" alt="DOI"></a>

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

Run Lynx:
```bash
lynx -d db.fa -o . example/*.fasta
```

# Lynx
Lynx: lightweight profiling of antibiotic resistance genes from short-read metagenomes

## Quick start
### Installation
```bash
conda install -c bioconda -c conda-forge xinehc::lynx
```

### Download database
```bash
wget -qN --show-progress https://zenodo.org/records/17536667/files/db.fa
```

### Download example
```bash
wget -qN --show-progress https://zenodo.org/records/15681353/files/example.tar.gz
tar -xvf example.tar.gz
```

### Run Lynx
```bash
lynx -d db.fa -o . example/*.fasta
```

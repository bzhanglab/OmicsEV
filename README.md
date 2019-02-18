## OmicsEV
A tool for large scale omics datasets evaluation

[<img src="https://github.com/bzhanglab/OmicsEV/blob/gh-pages/data/OmicsEV_overview.png" width=500 class="center">](https://bzhanglab.github.io/OmicsEV/)

## Installation

#### Method 1:

Install OmicsEV into the installed R library:
``` r
# Install the development version from GitHub:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
BiocManager::install("theislab/kBET")
BiocManager::install("wenbostar/metaX")
BiocManager::install("bzhanglab/NetSAM")
BiocManager::install("AnalytixWare/ShinySky")
BiocManager::install("bzhanglab/OmicsEV")
```
#### Method 2:
Install OmicsEV using Docker:
```sh
git clone -b docker https://github.com/bzhanglab/OmicsEV/
docker build -t omicsev OmicsEV
```

## Usage

Please follow the instruction in website *[https://bzhanglab.github.io/OmicsEV/](https://bzhanglab.github.io/OmicsEV/)*

## Example data

Please download an example dataset here: [RNA_seq_6_datasets.tar.gz](https://github.com/bzhanglab/OmicsEV/raw/gh-pages/data/RNA_seq_6_datasets.tar.gz)

## Example report

The HTML-based report for the example data: [OmicsEV report](https://bzhanglab.github.io/OmicsEV/data/example_report.html)

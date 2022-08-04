## OmicsEV


A tool for large scale omics data tables evaluation

[<img src="https://github.com/bzhanglab/OmicsEV/blob/gh-pages/data/OmicsEV_overview.png" width=500 class="center">](https://bzhanglab.github.io/OmicsEV/)

## Installation

#### Method 1 (recommended):
Install OmicsEV using Docker:
```sh
docker pull proteomics/omicsev
```

Use OmicsEV in docker:
```sh
# change the path your_data_path to a real path
docker run -it -v /your_data_path/:/opt/ -u $(id -u):$(id -g) proteomics/omicsev
```

#### Method 2:

Follow the instruction at https://github.com/bzhanglab/OmicsEV/wiki/OmicsEV-package-installation to install OmicsEV.


## Usage

Please follow the instruction in website *[https://bzhanglab.github.io/OmicsEV/](https://bzhanglab.github.io/OmicsEV/)*

## Examples

### Example 1: evaluate RNA-Seq data tables generated using different normalization methods

The RNA-Seq data is from TCGA-BRCA project. A total of six different data tables were generated using different normalization methods. A proteomics data table is available and it was generated from the same samples. Below is the R code to run the evaluation using OmicsEV. 

```r
library(OmicsEV)
run_omics_evaluation(data_dir = "datasets/",
                     sample_list = "sample_list.tsv",
                     x2 = "protein.tsv",
                     x2_label = "Protein",
                     cpu=0,
                     use_existing_data=TRUE,
                     data_type="gene",
                     class_for_ml="sample_ml.tsv")
```
Please download the input files for above code at [RNA_seq_6_datasets.tar.gz](https://github.com/bzhanglab/OmicsEV/raw/gh-pages/data/RNA_seq_6_datasets.tar.gz).

The HTML report generated using above code is available at [OmicsEV report](https://bzhanglab.github.io/OmicsEV/data/rna-seq-example_report.html).

This example takes about 2 hours and 40 minutes on a Linux system with 64 CPUs and 256G memory.


### Example 2: evaluate proteomics data tables generated using different pipelines

The proteomics data is from CPTAC Breast project. A total of three different data tables were generated using different pipelines. An RNA-Seq data table is available and it was generated from the same samples. Below is the R code to run the evaluation using OmicsEV. 

```r
library(OmicsEV)
run_omics_evaluation(data_dir = "datasets_75/",
                     sample_list = "sample_list_v2.tsv",
                     x2 = "rna.tsv",
                     cpu=0,
                     use_existing_data=TRUE,
                     data_type="gene",
                     class_for_ml="sample_ml.tsv")
```
Please download the input files for above code at [proteomics_3_datasets.tar.gz](https://github.com/bzhanglab/OmicsEV/raw/gh-pages/data/proteomics_3_datasets.tar.gz).

The HTML report generated using above code is available at [OmicsEV report](https://bzhanglab.github.io/OmicsEV/data/proteomics-example_report.html).

This example takes about one hour on a Linux system with 64 CPUs and 256G memory.





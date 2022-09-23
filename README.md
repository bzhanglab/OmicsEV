## OmicsEV


A tool for large scale omics data tables evaluation

[<img src="https://github.com/bzhanglab/OmicsEV/blob/gh-pages/data/OmicsEV_overview.png" width=500 class="center">](https://bzhanglab.github.io/OmicsEV/)

## Installation

Install OmicsEV using Docker:
```sh
docker pull proteomics/omicsev
```

If [Docker Engine](https://docs.docker.com/engine/) is not installed, please first install Docker Engine following the instruction at https://docs.docker.com/engine/install/. If [Docker Engine](https://docs.docker.com/engine/) is already installed, the OmicsEV docker can be installed using the above command line.

Use OmicsEV in docker:
```sh
# change the path your_data_path to a real path
docker run -it -v /your_data_path/:/opt/ -u $(id -u):$(id -g) proteomics/omicsev
# then lauch R 
R
```

Please put all the input data files for OmicsEV under a folder (for example: /your_data_path/, this can be any folder with write permission) and use parameter -v to map this folder to the Docker container directory "/opt/" (-v /your_data_path/:/opt/, don't change /opt/ part) so that all the input data files can be accessed inside OmicsEV docker. After lauching R in OmicsEV docker using above code, users can then use the OmicsEV functions to perform analysis. A few examples can be found below.

It requires a basic understanding of docker to use OmicsEV inside docker: https://www.docker.com/get-started/.


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
Please download the input files for above code at [RNA_seq_6_datasets.tar.gz](https://github.com/bzhanglab/OmicsEV/raw/gh-pages/data/RNA_seq_6_datasets.tar.gz). It contains the following files:

```shell
├── datasets
│   ├── d1.tsv
│   ├── d2.tsv
│   ├── d3.tsv
│   ├── d4.tsv
│   ├── d5.tsv
│   └── d6.tsv
├── protein.tsv
├── run_OmicsEV.R
├── sample_list.tsv
└── sample_ml.tsv
```

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
Please download the input files for above code at [proteomics_3_datasets.tar.gz](https://github.com/bzhanglab/OmicsEV/raw/gh-pages/data/proteomics_3_datasets.tar.gz). It contains the following files:

```
├── datasets_75
│   ├── CDAP.tsv
│   ├── MQ_ratio.tsv
│   └── paper.tsv
├── rna.tsv
├── run_OmicsEV.R
├── sample_list_v2.tsv
└── sample_ml.tsv
```

The HTML report generated using above code is available at [OmicsEV report](https://bzhanglab.github.io/OmicsEV/data/proteomics-example_report.html).

This example takes about one hour on a Linux system with 64 CPUs and 256G memory.


### Example 3: evaluate the data quality of single data table

The proteomics data is from CPTAC Breast project. A single data table was generated using one pipeline. An RNA-Seq data table is available and it was generated from the same samples. Below is the R code to run the evaluation using OmicsEV. 

```r
library(OmicsEV)
run_omics_evaluation(data_dir = "datasets_75/",
                     sample_list = "sample_list_v2.tsv",
                     x2 = "rna.tsv",
                     cpu=0,
                     use_existing_data=TRUE,
                     data_type="gene",
                     class_for_ml=c("LumA","LumB"))
```
Please download the input files for above code at [proteomics_1_dataset.tar.gz](https://github.com/bzhanglab/OmicsEV/raw/gh-pages/data/proteomics_1_dataset.tar.gz). It contains the following files:
```
├── datasets_75
│   └── paper.tsv
├── rna.tsv
├── run_OmicsEV.R
├── sample_list_v2.tsv
└── sample_ml.tsv
```

The HTML report generated using above code is available at [OmicsEV report](https://bzhanglab.github.io/OmicsEV/data/proteomics-one_dataset-example_report.html).

This example takes about one hour on a Linux system with 64 CPUs and 256G memory.



## Applications in publications
1. Cao, L., et al. [**Proteogenomic characterization of pancreatic ductal adenocarcinoma**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8654574/). *Cell* 2021;184(19):5031-5052 e5026.
2. Dou, Y., et al. [**Proteogenomic Characterization of Endometrial Carcinoma**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7233456/). *Cell* 2020;180(4):729-748 e726.
3. Gao, Q., et al. [**Integrated Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma**](https://pubmed.ncbi.nlm.nih.gov/31585088/). *Cell* 2019;179(2):561-577 e522.
4. Huang, C., et al. [**Proteogenomic insights into the biology and treatment of HPV-negative head and neck squamous cell carcinoma**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7946781/). *Cancer Cell* 2021;39(3):361-379 e316.
5. Satpathy, S., et al. [**Microscaled proteogenomic methods for precision oncology**](https://www.nature.com/articles/s41467-020-14381-2). *Nat Commun* 2020;11(1):532.
6. Anurag, M., et al. [**Proteogenomic markers of chemotherapy resistance and response in triple negative breast cancer**](https://doi.org/10.1158/2159-8290.CD-22-0200). *Cancer Discovery* 2022; CD-22.

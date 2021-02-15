# Install the development version of OmicsEV from GitHub:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
BiocManager::install("theislab/kBET")
BiocManager::install("wenbostar/metaX")
BiocManager::install("bzhanglab/NetSAM")
#BiocManager::install("AnalytixWare/ShinySky")
BiocManager::install("doMC")
BiocManager::install("bzhanglab/OmicsEV")


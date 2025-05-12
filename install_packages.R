# Run once from R to install the packages used by these scripts.
install.packages(c("statmod"), repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c("edgeR", "limma"), ask=FALSE, update=FALSE)

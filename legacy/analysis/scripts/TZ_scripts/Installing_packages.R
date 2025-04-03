library(remotes)
library(devtools)
library(magrittr)
library(decontX)

remotes::install_version(package = 'janitor', version = package_version('2.2.0'))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

install.packages("magrittr")
devtools::install_github("immunogenomics/harmony")
#remotes::install_version(package = 'harmony', version = package_version('1.1.0'))

remotes::install_version(package = 'Seurat', version = package_version('4.3.0'))
BiocManager::install("decontX", force = TRUE)

remotes::install_version("SeuratObject", version = "4.1.3")

install.packages("symphony")
devtools::install_github('satijalab/seurat-data@v0.2.1')
devtools::install_github("satijalab/azimuth@v0.4.6")
BiocManager::install("scRepertoire")

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))


remotes::install_github("satijalab/azimuth")
remotes::install_github("satijalab/seurat-data")


#in case of scDblFinder problems:
remove.packages("Matrix")
remove.packages("irlba")
install.packages("Matrix",type="source", dependencies=T)
install.packages("irlba",type="source", dependencies=T)

#packages for cell type annotation

install.packages("easybio", repos = c("https://person-c.r-universe.dev", "https://cloud.r-project.org"))



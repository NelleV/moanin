library(devtools)

if(!(require(pkgdown))){
    install.packages("pkgdown")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocCheck")

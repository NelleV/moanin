library(devtools)

if(!(require(pkgdown))){
    install.packages("pkgdown")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocCheck")

if(!(require(timecoursedata))){
    devtools::install_github("NelleV/timecoursedata")
}


if(!requireNamespace("BiocStyle", quietly=TRUE)){
    BiocManager::install("BiocStyle")
}


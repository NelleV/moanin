if(!require("devtools")){
    install.packages("devtools")
}

if(!require("testthat")){
    install.packages("testthat")
}

if(!require(lattice)){
    install.packages("lattice")
}


if(!require(knitr)){
    install.packages("knitr")
}

if(!require(rmarkdown)){
    install.packages("rmarkdown")
}

if(!require(roxygen2)){
    install.packages("roxygen2")
}


if(!(require(pkgdown))){
    install.packages("pkgdown")
}

if (!require("BiocManager")){
    install.packages("BiocManager")
}

if(!(require(BiocCheck))){
    BiocManager::install("BiocCheck")
}

if(!(require(viridis))){
    install.packages("viridis")
}

if(!(require(reshape2))){
    install.packages("reshape2")
}

if(!require(testthat)){
    install.packages("testthat")
}


if(!require(roxygen2)){
    install.packages("roxygen2")
}

if(!require(ClusterR)){
    install.packages("ClusterR")
}

if(!require(splines)){
    install.packages("splines")
}


if(!require(limma)){
    BiocManager::install("limma")
}

if(!require(edgeR)){
    BiocManager::install("edgeR")
}

if(!require(edge)){
    BiocManager::install("edge")
}

if(!require(ClusterR)){
    install.packages("ClusterR")
}

if(!require(covr)){
    install.packages("covr")
}

if(!require(kableExtra)){
    install.packages("kableExtra")
}

if(!require(shiny)){
    install.packages("shiny")
}

if(!require(DT)){
    install.packages("DT")
}

if(!require(topGO)){
    BiocManager::install("topGO")
}

if(!require(MASS)){
    BiocManager::install("MASS")
}


if(!require("DelayedArray")){
    BiocManager::install("DelayedArray")
}

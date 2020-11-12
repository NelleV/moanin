#' Provides set of basis functions on either side of a time point, allowing for a discontinuity in the fitted functions
#' @param timepoints vector of numeric timepoints for which the splines basis will be evaluated
#' @param discon_point a single numeric value that represents where the discontinuity should be
#' @param knots passed to \code{ns} or \code{bs}. If not NULL, should give knots on either side of \code{discon_point} as single vector -- they will be separated in the call to \code{discon_point}
#' @param intercept Whether to include an intercept (vector of all 1s) for each side of the discontinuity. Note this is different than the argument \code{intercept} of either \code{bs} or \code{ns}, which is set to \code{FALSE}. 
#' @param type either "ns" or "bs" indicating which splines basis function should be used.
#' @param degree passed to \code{bs} (if applicable)
#' @param dfPre the df for the basis functions defined before the discontinuity point
#' @param dfPost the df for the basis functions defined after the discontinuity point
#' @examples
#' x<-seq(0,10,length=100)
#' basis<-discont_basis(x,discont_point=3, dfPre=3, dfPost=4, intercept=TRUE)
#' # Plot of the basis functions
#' par(mfrow=c(3,3))
#' for(i in 1:ncol(basis)){
#'    plot(x,basis[,i],type="l")
#'    abline(v=3,lty=2)
#' }
#' # Use it in a moanin_model object instead of ns/bs:
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData, meta=testMeta,
#'     spline_formula=~Group:discont_basis(Timepoint,dfPre=3,
#'         dfPost=3,discont=20,intercept=TRUE)+0,
#'     df=6)

discont_basis<-function(timepoints, discont_point, 
    knots=NULL, dfPre=NULL, dfPost=dfPre,degree=3, 
        intercept=TRUE,type=c("ns","bs")){

        # Useful for definition of intercept: https://cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf
        # Basically not the flat line intercept, but will create perfect multi-coliniearity if include it because for every observation will sum to 1 across the basis functions.
        
        
        #The df in command bs actually means the number of columns of the design matrix returned by bs.So if the intercept is not included in the design matrix (which is the default), then the df in command bs is equal to the real df minus 1....If intercept = TRUE, then we need m=df-2 knots, otherwise we need m=df-1 knots. http://publish.illinois.edu/liangf/files/2016/05/Note_splines.pdf

    if(!is.null(dfPre) & is.null(dfPost)) dfPost<-dfPre
    if(is.null(dfPre) & !is.null(dfPost)) dfPre<-dfPost        
    if(is.null(dfPre)& is.null(dfPost) & is.null(knots)) 
        stop("Must provide either dfPre and dfPost or knots")
    type<-match.arg(type)
    prePoints<-which(timepoints<=discont_point)                     
    postPoints<-which(timepoints>discont_point)    
    if(!is.null(knots)){
        knotsPre<-knots[knots<=discont_point]                 
        knotsPost<-knots[knots>discon_point]
    }
    else{
        knotsPre<-knotsPost<-NULL
    }
    bsfun<-switch(type,
        "ns"=function(...){splines::ns(...)},
        "bs"=function(...){splines::bs(degree=degree,...)}
        )  
    preBasis<-bsfun(x=timepoints[prePoints],
          knots=knotsPre,df=dfPre,
          intercept=FALSE)
    postBasis<-bsfun(x=timepoints[postPoints],
        knots=knotsPost,df=dfPost,
        intercept=FALSE)
                         
    if(intercept){
        preBasis<-cbind(rep(1,nrow(preBasis)),preBasis)
        postBasis<-cbind(rep(1,nrow(postBasis)),postBasis)
    }
    combinedBasis<-matrix(0,
        nrow=length(timepoints),
        ncol=ncol(preBasis)+ncol(postBasis))     
    
    preCols<-1:ncol(preBasis)
    postCols<-(ncol(preBasis)+1):(ncol(combinedBasis))
    combinedBasis[prePoints,preCols]<-preBasis
    combinedBasis[postPoints,postCols]<- postBasis
    return(combinedBasis)   
                      
}


# ## Test EPICON Points:
# dataDir<-"/Users/epurdom/Documents/Research/EPICON/epiconWorkRepos/mRNASeq/scripts/results/Year1"
# pybasis<-read.csv(file.path(dataDir,"basis/leaf_basis_preflowering.txt"))
# barcodes<-pybasis[,1]
# meta<-read.table(file.path(dataDir,"data/leaf_meta.tsv"))
# meta<-meta[match(barcodes,meta$Barcode),]
# fac<-factor(meta$Genotype):factor(meta$Condition)
#
# pybasis<-data.matrix(pybasis[,-1])
# pyord<-do.call("order",as.data.frame(t(pybasis)))
# pybasis<-pybasis[,pyord]
#
#
# #bs says degrees of freedom is too small -- makes it 4
# basis<-create_basis(meta$Time.Point,discont_point=8.5,
#     dfPre=3, dfPost=3,degree=3, intercept=TRUE,type="bs")
# rbasis<-model.matrix(~Genotype:Condition:basis+0,data=meta)
# rbasis<-unname(rbasis)
# rord<-do.call("order",as.data.frame(t(rbasis)))
# rbasis<-rbasis[,rord]
#
# # basisNoInt<-create_basis(meta$Time.Point,discont_point=8.5,
# #     dfPre=3, dfPost=3,degree=3, intercept=FALSE,type="bs")
# # prepostFac<-factor(meta$Time.Point<=8.5)
# # rbasisNI<-model.matrix(~Genotype:Condition:basisNoInt+Genotype:Condition:prepostFac+0,data=meta)
# # rordNI<-do.call("order",as.data.frame(t(rbasisNI)))
# # rbasisNI<-rbasisNI[,rordNI]
# #rbasisNI<-unname(rbasisNI)
# dim(rbasis)
# dim(pybasis)
# dim(rbasisNI)
#
#
# # simTimepoints<-rep(seq(1,17,length=100),2)
# # simCondition<-sort(gl(n=2,k=2,length=100))
# # simBasis<-model.matrix(~simCondition:create_basis(simTimepoints,discont_point=8.5,
# #     dfPre=3, dfPost=3,degree=3, intercept=TRUE,type="bs")
#
#
# ##########
# ##Get same order
# ##########
# #sort by factor
#
#
# plotPerLevel<-function(level,mat,basisCol,add=TRUE,col="black"){
#    wh<-which(fac==level)
#    ord<-order(meta$Time.Point[wh])
#    if(add)
#        lines(meta$Time.Point[wh][ord],
#            mat[wh,basisCol][ord],type="l",lty=1,col=col)
#    else plot(meta$Time.Point[wh][ord],
#         mat[wh,basisCol][ord],type="l",lty=1,col=col,
#         xlab="Timepoint",ylab="Basis",ylim=c(0,1))
#
# }
# levs<-levels(fac)
# dev.set(2)
# par(mfrow=c(6,6))
# for(i in 1:ncol(rbasis)){
#    for(ll in 1:length(levs)){
#      plotPerLevel(level=levs[ll],mat=rbasis,basisCol=i,add=(ll!=1), col=palette()[ll])
#    }
#    if(i==3) title(main="R Basis")
#
#    if(i==ncol(rbasis)){
#        plot.new()
#        legend("topleft",levs,fill=palette()[1:4])
#    }
#
# }
# dev.set(2)
# par(mfrow=c(6,7))
# for(i in 1:ncol(rbasisNI)){
#    for(ll in 1:length(levs)){
#      plotPerLevel(level=levs[ll],mat=rbasisNI,basisCol=i,add=(ll!=1), col=palette()[ll])
#    }
#    if(i==3) title(main="R Basis, different Int")
#
#    if(i==ncol(rbasis)){
#        plot.new()
#        legend("topleft",levs,fill=palette()[1:4])
#    }
#
# }
# dev.set(3)
# par(mfrow=c(6,6))
# for(i in 1:ncol(pybasis)){
#    for(ll in 1:length(levs)){
#      plotPerLevel(level=levs[ll],mat=pybasis,basisCol=i,add=(ll!=1), col=palette()[ll])
#
#    }
#    if(i==3) title(main="Python Basis")
#    if(i==ncol(pybasis)){
#        plot.new()
#        legend("topleft",levs,fill=palette()[1:4])
#    }
#
# }




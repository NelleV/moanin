# Returns the basis matrix
#
# Parameters
# ----------
# timepoints : ndarray (t, ) of timepoints
#
# conditions : ndarray (t, ), optional, default: None
#     ndarray containing the conditions information.
#
# knots : ndarray, optional, default: None
#     The interior knots to use for the spline. If unspecified, then equally
#     spaced quantiles of the input data are used. You must specify at least
#     one of ``df`` and ``knots``
#
# df : int, default 4
#     Number of degrees of freedom of the spline. This corresponds to the
#     size of the basis
#
# degree : int, default: 3
#     The degree of the splines to use
#
# Returns
# -------
# basis : ndarray (t, k)
# """

create_basis<-function(timepoints, discont_point, 
    knots=NULL, dfPre=NULL, dfPost=dfPre,degree=3, orthonormalize=FALSE,
    intercept=FALSE,type=c("ns","bs")){
    if(!is.null(dfPre) & is.null(dfPost)) dfPost<-dfPre
    if(is.null(dfPre) & !is.null(dfPost)) dfPre<-dfPost        
    if(is.null(dfPre)& is.null(dfPost) & is.null(knots)) 
        stop("Must provide either dfPre and dfPost or knots")
    type<-match.arg(type)
    prePoints<-which(timepoints<=discont_point)                     
    postPoints<-which(timepoints>discont_point)                     
    bsfun<-switch(type,
        "ns"=function(...){splines::ns(...)},
        "bs"=function(...){splines::bs(degree=degree,...)}
        )  
    preBasis<-bsfun(x=timepoints[prePoints],
          knots=knots,df=dfPre,
          intercept=intercept)
    postBasis<-bsfun(x=timepoints[postPoints],
        knots=knots,df=dfPost,
        intercept=intercept)
                         
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
# pybasis<-data.matrix(pybasis[,-1])
# meta<-read.table(file.path(dataDir,"data/leaf_meta.tsv"))
# meta<-meta[match(barcodes,meta$Barcode),]
# #bs says degrees of freedom is too small -- makes it 4
# basis<-create_basis(meta$Time.Point,discont_point=8.5,
#     dfPre=4, dfPost=4,degree=3, intercept=TRUE,type="bs")
# rbasis<-model.matrix(~Genotype:Condition:basis+0,data=meta)
# rbasis<-unname(rbasis)
# dim(rbasis)
# dim(pybasis) #First column is barcode
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
# fac<-factor(meta$Genotype):factor(meta$Condition)
#
# pyord<-do.call("order",as.data.frame(t(pybasis)))
# rord<-do.call("order",as.data.frame(t(rbasis)))
# head(pybasis[,pyord])
# head(rbasis[,rord])
#
#
# pybasis<-pybasis[,pyord]
# rbasis<-rbasis[,rord]
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




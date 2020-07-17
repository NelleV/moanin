

setGeneric("log_transform", 
           function(object) { standardGeneric("log_transform")})
setGeneric("get_log_data", 
           function(object) { standardGeneric("get_log_data")})
setGeneric("degrees_of_freedom", 
    function(object) { standardGeneric("degrees_of_freedom")})
setGeneric("basis_matrix", 
    function(object) { standardGeneric("basis_matrix")})
setGeneric("spline_formula", 
    function(object) { standardGeneric("spline_formula")})
setGeneric("group_variable_name", 
    function(object) { standardGeneric("group_variable_name")})
setGeneric("group_variable_name<-", 
    function(object, value) { standardGeneric("group_variable_name<-")})
setGeneric("group_variable", 
    function(object) { standardGeneric("group_variable")})
setGeneric("group_variable<-", 
    function(object, value) { standardGeneric("group_variable<-")})
setGeneric("time_variable_name", 
    function(object) { standardGeneric("time_variable_name")})
setGeneric("time_variable_name<-", 
    function(object,  value) { standardGeneric("time_variable_name<-")})
setGeneric("time_variable", 
    function(object) { standardGeneric("time_variable")})
setGeneric("time_variable<-", 
    function(object,  value) { standardGeneric("time_variable<-")})
setGeneric("time_by_group_variable", 
    function(object) { standardGeneric("time_by_group_variable")})

#' @name Moanin-methods
#' @title Helper methods for the Moanin class
#'
#' @description This is a collection of helper methods for the Moanin
#'   class.
#' @inheritParams DE_timecourse
#' @param value replacement value
#' @param x \code{Moanin} object
#' @param ... arguments passed to subsetting
#' @examples 
#' # Load some data
#' data(exampleData)
#' moanin = create_moanin_model(data=testData,meta=testMeta)
#' group_variable_name(moanin)
#' time_variable_name(moanin)
#' @return \code{group_variable_name} and \code{time_variable_name} return the 
#' name of the column containing the variable. \code{group_variable} and 
#' \code{time_variable} return the  actual variable. 
#' @export
#' @aliases group_variable_name group_variable_name,Moanin-method
setMethod("group_variable_name","Moanin",function(object){
    object@group_variable_name
})
#' @rdname Moanin-methods
#' @aliases group_variable_name<- group_variable_name<-,Moanin-method
#' @export
setReplaceMethod("group_variable_name","Moanin",function(object,value){
    object@group_variable_name<-value
    return(object)
    
})
#' @rdname Moanin-methods
#' @aliases time_by_group_variable time_by_group_variable,Moanin-method
#' @export
setMethod("time_by_group_variable","Moanin",function(object){
    colData(object)$WeeklyGroup
})
#' @rdname Moanin-methods
#' @aliases group_variable group_variable,Moanin-method
#' @export
setMethod("group_variable","Moanin",function(object){
    gpVar<-group_variable_name(object)
    colData(object)[,gpVar]
})
#' @rdname Moanin-methods
#' @aliases group_variable<- group_variable<-,Moanin-method
#' @export
setReplaceMethod("group_variable","Moanin",function(object,value){
    gpVar<-group_variable_name(object)
    colData(object)[,gpVar]<-value
    return(object)
    
})
#' @rdname Moanin-methods
#' @aliases time_variable_name time_variable_name,Moanin-method
#' @export
setMethod("time_variable_name","Moanin",function(object){
    object@time_variable_name
})
#' @rdname Moanin-methods
#' @aliases time_variable_name<- time_variable_name<-,Moanin-method
#' @export
setReplaceMethod("time_variable_name","Moanin",function(object,value){
    object@time_variable_name<-value
    return(object)
    
})
#' @rdname Moanin-methods
#' @aliases time_variable time_variable,Moanin-method
#' @export
setMethod("time_variable","Moanin",function(object){
    tpVar<-time_variable_name(object)
    colData(object)[,tpVar]
})
#' @rdname Moanin-methods
#' @aliases time_variable<- time_variable_name<-,Moanin-method
#' @export
setReplaceMethod("time_variable","Moanin",function(object,value){
    tpVar<-time_variable_name(object)
    colData(object)[,tpVar]<-value
    return(object)
})
#' @rdname Moanin-methods
#' @aliases degrees_of_freedom degrees_of_freedom,Moanin-method
#' @export
setMethod("degrees_of_freedom","Moanin",function(object){
    object@degrees_of_freedom
})
#' @rdname Moanin-methods
#' @aliases basis_matrix basis_matrix,Moanin-method
#' @export
setMethod("basis_matrix","Moanin",function(object){
    object@basis_matrix
})
#' @rdname Moanin-methods
#' @aliases spline_formula spline_formula,Moanin-method
#' @export
setMethod("spline_formula","Moanin",function(object){
    object@spline_formula
})
#' @rdname Moanin-methods
#' @aliases show show,Moanin-method
#' @export
setMethod("show","Moanin",function(object){
        N<-ncol(object)
        cat("Moanin object on",N,
            "samples containing the following information:\n")
        cat(paste0("Group variable given by '",group_variable_name(object),
                   "' with the following levels:\n"))
        print(summary(group_variable(object)))
        cat(paste0("Time variable given by '",time_variable_name(object),"'\n"))
        cat("Basis matrix with",ncol(basis_matrix(object)),
            "basis_matrix functions\n")
        if(!is.null(spline_formula(object))){
            cat("Basis matrix was constructed with", 
                "the following spline_formula\n")
            #get rid of extra spaces:
            form<-gsub("\\s{2,}","",deparse(spline_formula(object))) 
            cat(paste(form,collapse="",sep=""),"\n")
        }
        else{
            if(is.null(degrees_of_freedom(object)))
                cat("Basis matrix was provided by user,",
                "spline_formula and degrees_of_freedom=NULL\n")
            else
                cat("Basis matrix and degrees of freedom provided by user,",
                    "equal to", degrees_of_freedom(object),"\n")
        }
        cat("\nInformation about the data (a SummarizedExperiment object):\n")
        show(as(object,"SummarizedExperiment"))
})


### For subsetting
#' @details Note that when subsetting the data, the dendrogram information and
#'   the co-clustering matrix are lost.
#' @aliases [,Moanin,ANY,ANY,ANY-method [,Moanin,ANY,character,ANY-method
#' @param i,j A vector of logical or integer subscripts, indicating the rows and
#'   columns to be subsetted for \code{i} and \code{j}, respectively.
#' @param drop A logical scalar that is ignored.
#' @rdname Moanin-methods
#' @export
setMethod(
    f = "[",
    signature = c("Moanin", "ANY", "character"),
    definition = function(x, i, j, ..., drop=TRUE) {
        j<-match(j, colnames(x))
        callGeneric()
        
    }
)
#' @rdname Moanin-methods
#' @export
setMethod(
    f = "[",
    signature = c("Moanin", "ANY", "logical"),
    definition = function(x, i, j, ..., drop=TRUE) {
        j<-which(j)
        callGeneric()
    }
)
#' @rdname Moanin-methods
#' @export
setMethod(
    f = "[",
    signature = c("Moanin", "ANY", "numeric"),
    definition = function(x, i, j, ..., drop=TRUE) {

        out<- new("Moanin",
                #have to explicitly give the inherintence... not great:
                as(selectMethod("[",
                    c("SummarizedExperiment","ANY","numeric"))(x,i,j),
                   "SummarizedExperiment"),
                basis_matrix=basis_matrix(x)[j, , drop=FALSE],
                spline_formula=spline_formula(x),
                degrees_of_freedom=degrees_of_freedom(x),
                group_variable_name=group_variable_name(x),
                time_variable_name=time_variable_name(x)
        )
        return(out)
    }
)

#' return log data if appropriate
#' @rdname Moanin-methods
#' @aliases log_transform log_transform,Moanin-method
#' @export
setMethod("log_transform","Moanin",function(object){
    object@log_transform
})

#' return log data if appropriate
#' @rdname Moanin-methods
#' @aliases get_log_data get_log_data,Moanin-method
#' @export
setMethod("get_log_data","Moanin",function(object){
    if(log_transform(object)) y<-log(assay(object)+1)
    else y<-assay(object)
    return(data.matrix(y))
})
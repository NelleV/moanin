setGeneric("df", 
           function(object) { standardGeneric("df")})
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

#' @name Moanin-methods
#' @title Helper methods for the Moanin class
#'
#' @description This is a collection of helper methods for the Moanin
#'   class.
#' @examples 
#' # Load some data
#' data(exampleData)
#' moanin = create_moanin_model(data=testData,meta=testMeta)
#' group_variable_name(moanin)
#' time_variable_name(moanin)

setMethod("group_variable_name","Moanin",function(object){
    object@group_variable_name
})
#' @rdname Moanin-methods
setReplaceMethod("group_variable_name","Moanin",function(object,value){
    object@group_variable_name<-value
    return(object)
    
})

#' @rdname Moanin-methods
setMethod("group_variable","Moanin",function(object){
    gpVar<-group_variable_name(object)
    colData(object)[,gpVar]
})
setReplaceMethod("group_variable","Moanin",function(object,value){
    gpVar<-group_variable_name(object)
    colData(object)[,gpVar]<-value
    return(object)
    
})
#' @rdname Moanin-methods
setMethod("time_variable_name","Moanin",function(object){
    object@time_variable_name
})
#' @rdname Moanin-methods
setReplaceMethod("time_variable_name","Moanin",function(object,value){
    object@time_variable_name<-value
    return(object)
    
})
#' @rdname Moanin-methods
setMethod("time_variable","Moanin",function(object){
    tpVar<-time_variable_name(object)
    colData(object)[,tpVar]
})
#' @rdname Moanin-methods
setReplaceMethod("time_variable","Moanin",function(object,value){
    tpVar<-time_variable_name(object)
    colData(object)[,tpVar]<-value
    return(object)
})
#' @rdname Moanin-methods
setMethod("df","Moanin",function(object){
    object@degrees_of_freedom
})
#' @rdname Moanin-methods
setMethod("basis_matrix","Moanin",function(object){
    object@basis_matrix
})
#' @rdname Moanin-methods
setMethod("spline_formula","Moanin",function(object){
    object@spline_formula
})
#' @rdname Moanin-methods
setMethod("show","Moanin",function(object){
        N<-ncol(object)
        cat("Moanin object on",N,"samples containing the following information:\n")
        cat(paste0("Group variable given by '",group_variable_name(object),"' with the following levels:\n"))
        print(summary(group_variable(object)))
        cat(paste0("Time variable given by '",time_variable_name(object),"'\n"))
        cat("Basis matrix with",ncol(basis_matrix(object)),"basis_matrix functions\n")
        if(!is.null(spline_formula(object))){
            cat("Basis matrix was constructed with the following spline_formula\n")
            form<-gsub("\\s{2,}","",deparse(spline_formula(object))) #get rid of extra spaces
            cat(paste(form,collapse="",sep=""),"\n")
        }
        else{
            if(is.null(df(object)))
                cat("Basis matrix was provided by user, spline_formula and degrees_of_freedom=NULL\n")
            else
                cat("Basis matrix and degrees of freedom provided by user, equal to",df(object),"\n")
        }
    
})
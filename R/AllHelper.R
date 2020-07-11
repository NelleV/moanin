setGeneric("group_variable_name", 
           function(object,  ...) { standardGeneric("group_variable_name")})
setGeneric("group_variable_name<-", 
           function(object, value) { standardGeneric("group_variable_name")})
setGeneric("group_variable", 
           function(object,  ...) { standardGeneric("group_variable")})
setGeneric("group_variable<-", 
           function(object, value) { standardGeneric("group_variable")})
setGeneric("time_variable_name", 
           function(object,  ...) { standardGeneric("time_variable_name")})
setGeneric("time_variable_name<-", 
           function(object,  value) { standardGeneric("time_variable_name")})
setGeneric("time_variable", 
           function(object,  ...) { standardGeneric("time_variable")})
setGeneric("time_variable<-", 
           function(object,  value) { standardGeneric("time_variable")})

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

setMethod("group_variable_name","Moanin",function(object,...){
    object@group_variable
})
#' @rdname Moanin-methods
setReplaceMethod("group_variable_name","Moanin",function(object,value){
    object@group_variable<-value
    return(object)
    
})

#' @rdname Moanin-methods
setMethod("group_variable","Moanin",function(object,...){
    gpVar<-group_variable_name(object)
    colData(object)[,gpVar]
})
setReplaceMethod("group_variable","Moanin",function(object,value){
    gpVar<-group_variable_name(object)
    colData(object)[,gpVar]<-value
    return(object)
    
})
#' @rdname Moanin-methods
setMethod("time_variable_name","Moanin",function(object,...){
    object@time_variable
})
#' @rdname Moanin-methods
setReplaceMethod("time_variable_name","Moanin",function(object,value){
    object@time_variable<-value
    return(object)
    
})
#' @rdname Moanin-methods
setMethod("time_variable","Moanin",function(object,...){
    tpVar<-time_variable_name(object)
    colData(object)[,tpVar]
})
#' @rdname Moanin-methods
setReplaceMethod("time_variable","Moanin",function(object,value){
    tpVar<-time_variable_name(object)
    colData(object)[,tpVar]<-value
    return(object)
})
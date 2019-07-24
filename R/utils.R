
# This should override the Biobase function
rowMax = function(data){
    check_is_2d(data)
    return(as.vector(apply(data, 1, max)))
}

rowMin = function(data){
    check_is_2d(data)
    return(as.vector(apply(data, 1, min)))
}


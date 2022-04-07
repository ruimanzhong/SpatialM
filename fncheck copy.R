fnCheckInputsMelding = function(depoint, dearea, dppoint, dparea) {
  
  if (is.null(dearea) == T && is.null(depoint) == T) stop("'Valid estimation data input required'")
  if (is.null(dparea) == T && is.null(dppoint) == T) stop("'Valid preidction data input required'")
  
  if(sum(c(class(depoint)[[1]],class(dppoint)[[1]],
           class(dparea)[[1]],class(dearea)[[1]]) %in% c("sf", "NULL")) != 4)
    stop("'All input data should be 'sf' obj")
  
  
  if(is.null(depoint) == F && sum(c("pvalue", "geometry") %in% colnames(depoint)) != 2 ) 
    stop("'depoint' must have 'geometry','pvalue' as column names")
  if(is.null(dearea) == F && sum(c("avalue", "geometry") %in% colnames(dearea)) != 2 ) 
    stop("'dearea' must have 'avalue','geometry', as column names")
  
  if(st_crs(depoint) != st_crs(dearea) && sum(c(is.null(depoint),is.null(dearea))) == 0){
    stop('all the input data must have the same crs, use st_crs() to check your data')
  }
   if(is.null(dparea) == F  && sum(c(st_crs(dparea) == st_crs(depoint),st_crs(dparea) == st_crs(dearea))) == 0){
     stop('all the input data must have the same crs, use st_crs() to check your data')
   }
  if(is.null(dppoint) == F && sum(c(st_crs(dppoint) == st_crs(depoint),st_crs(dppoint) == st_crs(dearea))) == 0){
    stop('all the input data must have the same crs, use st_crs() to check your data')
  }
}

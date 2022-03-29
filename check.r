check_meld = function(area.sf, p.sf, area.pre, proN)
{
  packages = c('sf','INLA')
  
  package.check = lapply(packages, FUN =function(x){
    if(!require(x, character.only = T))
      install.packages(x)
    if(!(x %in% ( .packages()  ) ) )
      library(x ,character.only = T)})  
  
  if (is.null(area.pre) == T ||  class(area.pre)[[1]] != 'sf') 
    stop("'Valid predict area input required, should be sf class'")
  if (is.null(area.sf) == T && is.null(p.sf) == T) stop("'Valid data input required'")
  if (is.null(area.sf) == F &&  class(area.sf)[[1]] != 'sf') stop("'area.sf' must be sf obj")
  if (is.null(p.sf) == F &&  class(p.sf) != 'sf') stop("'p.sf' must be sf obj")
  
  if(is.null(p.sf) == F && sum(c("value.p", "geometry") %in% colnames(p.sf)) != 2 ) 
    stop("'p.df' must have 'geometry','value.p' as column names")
  if(is.null(area.sf) == F && sum(c("value.a", "geometry") %in% colnames(area.sf)) != 2 ) 
    stop("'area.sf' must have 'value.a','geometry', as column names")
  
  if(is.null(area.sf) == F && st_crs(area.sf) != st_crs(area.pre) )
    stop("The crs of 'area.sf' and 'area.pre' is not the same")
  if(is.null(p.sf) == F && st_crs(p.sf) != st_crs(area.pre))
    stop("The crs of 'p.sf' and 'area.pre' is not the same'")

}

check_mesh = function(bd.sf) {
  packages = c('sf','INLA')
  
  package.check = lapply(packages, FUN =function(x){
    if(!require(x ,character.only = T))
      install.packages(x)
    if(!(x %in% ( .packages()  ) ) )
      library(x ,character.only = T)})  
  
  if (is.null(bd.sf)== T ||  class(bd.sf)[[1]] != 'sf') 
    stop("'bd.sf' should be a sf obj to construct a mesh")
}

check_dwscaler = function(area.sf, p.sf, area.pre)
{
  packages = c('sf','INLA')
  
  package.check = lapply(packages, FUN =function(x){
    if(!require(x ,character.only = T))
      install.packages(x)
    if(!(x %in% ( .packages()  ) ) )
      library(x ,character.only = T)})  
  
  if (class(area.sf)[[1]] != 'sf') 
    stop("'area.sf' must be sf obj")
  if (class(p.sf)[[1]] != 'sf') 
    stop("'p.sf' must be sf obj")
  if (class(area.pre)[[1]] != 'sf') 
    stop("'area.pre' must be sf obj")
  
  if( sum(c("value.p", "geometry") %in% colnames(p.sf)) != 2 ) 
    stop("'p.df' must have 'geometry','value.p' as column names")
  if( sum(c("value.a", "geometry") %in% colnames(area.sf)) != 2 ) 
    stop("'area.sf' must have 'value.a','geometry', as column names")
}

check_projM = function(Aa, AP){
  
  message('The number of the NULL point and areal matrix: ')
  return(is.null(Aa) + is.null(Ap))
  
  
}

check_pre = function(bd.sf =  NULL,pre.sf = NULL, pre.p = NULL, 
                     cell.x = NULL, cell.y = NULL, proN = 4326) {
  
  if (sum(is.null(bd.sf,pre.sf,pre.p)) != 2)
    stop("'bd.sf' is for surface prediction, 'pre.sf' is for areal prediction, and 
         'pre.p' is for point prediction, input one and only one of them")
  if(is.null(bd.sf) == F && (is.null(cell.x) == T |is.null(cell.y == T) ))
    stop("'bd.sf' is for surface prediction, 'cell.x' and 'cell.y' must be numbers to construct grids")
  if (is.null(pre.sf) == F ||  class(pre.sf)[[1]] != 'sf') 
    stop("'Valid predict area input required, ‘pre.sf’ should be sf class'")
  if (is.null(pre.p) == F ||  class(pre.p)[[1]] != 'sf') 
    stop("'Valid predict point input required, ‘pre.p’ should be sf class'")
}

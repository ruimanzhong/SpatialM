#' Target Area Construction
#' 
#' @param pre.sf A sf obj of the boundary of the area target area
#' @param cell.x A number to generate cells in the target area
#' @param cell.y A number to generate cells in the target area 
#' @param proN A number that gives a way to project the location. By default is 4326, which is WG84
#' @return A sf obj of target area 

source('check.r')
target = function( bd.sf = NULL, pre.sf = NULL, cell.x = NULL, 
                   cell.y = NULL,  proN = 4326) {
  
  if(is.null(pre.sf) == F && is.null(bd.sf) == T) {
    bd.sf = st_as_sf(st_union(pre.sf))
  }

  bd.sf = st_transform(bd.sf, proN)
  bb = unname(attributes(st_geometry(bd.sf))$bbox)
  x = seq(bb[1] - 1, bb[3] + 1, length.out = cell.x)
  y = seq(bb[2] - 1, bb[4] + 1, length.out = cell.y)
  coop = expand.grid(x, y)
  coop_sf = sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = proN)
  
  target = coop_sf %>% 
    st_join(bd.sf, left = FALSE)
  
  return(target)
  
  
}
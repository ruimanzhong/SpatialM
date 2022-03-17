#' Target Area Construction
#' 
#' @param pre.sf A sf obj of the boundary of the area target area
#' @param cell.x A number to generate cells in the target area
#' @param cell.y A number to generate cells in the target area 
#' @param proN A number that gives a way to project the location. By default is 4326, which is WG84
#' @return A sf obj of target area 


area.pre = function(pre.sf, cell.x, cell.y, proN = 4326) {
  
  pre.sf = st_transform(pre.sf, proN)
  bb = unname(attributes(st_geometry(pre.sf))$bbox)
  x = seq(bb[1] - 1, bb[3] + 1, length.out = cell.x)
  y = seq(bb[2] - 1, bb[4] + 1, length.out = cell.y)
  coop = expand.grid(x, y)
  coop_sf = sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = proN)
  
  pre.sf = pre.sf%>%
  st_set_crs(proN)
  
  return(
  coop_sf %>% 
    st_join(pre.sf, left = FALSE)
  )
}

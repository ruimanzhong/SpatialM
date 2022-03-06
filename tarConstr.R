#**************construct Prediction Area**********

area.pre = function(pre.sf,cell.x,cell.y,proN){
  bb<- unname(attributes(st_geometry(pre.sf))$bbox)
  x <- seq(bb[1] - 1, bb[3] + 1, length.out = cell.x)
  y <- seq(bb[2] - 1, bb[4] + 1, length.out = cell.y)
  coop <- expand.grid(x, y)
  coop_sf = sf::st_as_sf(coop,coords = c('Var1','Var2'),crs = proN)
  
  pre.sf = pre.sf%>%
  st_set_crs(proN)
  
  coop_sf %>% 
    st_join(pre.sf,left = FALSE)
}

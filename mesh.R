source('check.r')
Mesh = function(area.sf,offset = NULL,cut.off = 0.1) {
  
  check_mesh(area.sf)
  area.sf$centroids = st_transform(area.sf, proN) %>% 
    st_centroid() 
  loc.sf = as.matrix(st_coordinates(area.sf$centroids[,1]))[,c(1,2)]
  bd = st_union(area.sf)
  bd.ss  = as(bd,Class = 'Spatial')
  
  max.edge = max(attributes(st_geometry(area.sf))$bbox[3] -attributes(st_geometry(area.sf))$bbox[1],
                 attributes(st_geometry(area.sf))$bbox[4] -attributes(st_geometry(area.sf))$bbox[2])/10
 return(
  mesh = inla.mesh.2d(
    loc = loc.sf,
   boundary = bd.ss,
   max.edge = max.edge, cutoff = cutoff,
   offset = offset)
   ) 
}

fnCreateMesh = function(depoint, dearea, bd) {
  
  de1 <- depoint
  de2 <- dearea
  de1ToF <- !is.null(depoint)
  de2ToF <- !is.null(dearea)
  
  location = NULL
    
  if(de1ToF){ 
    de1 = st_transform(de1,crs = 4326)
    location = as.matrix(st_coordinates(de1)[ , c(1,2)])
  }
  
  bdsp  = as(bd,Class = 'Spatial')
  range =  fnRange(bd)
  
  return(
    inla.mesh.2d(
      loc = location,
      boundary = bdsp,
      max.edge = c(range/10, range), cutoff = range/25
    )
  )
}

fnRange = function(bd){
  
  range = 0.33*max(attributes(st_geometry(bd))$bbox[3] -attributes(st_geometry(bd))$bbox[1],
                   attributes(st_geometry(bd))$bbox[4] -attributes(st_geometry(bd))$bbox[2])
  return(range)
}
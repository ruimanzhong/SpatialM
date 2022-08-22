punifsample <- function(n,rdata, crsproj){
  dp <- rpoint(n)
  #plot(dp)
  # Extract data at points
  s <- as.data.frame(cbind(x = dp$x,y = dp$y, value = raster::extract(rdata, cbind(dp$x, dp$y))))
  depoint = s %>% st_as_sf(coords = c("x", "y"), dim = "XY") %>%
    st_set_crs(crsproj) %>% st_cast("MULTIPOINT")
  
  return(depoint)
}
pimsample <- function(n,rdata, crsproj){
  im_r <- as.im(rdata)
  # sample from im_r based on value of im_r
  dp <- rpoint(n,im_r)
  #plot(dp)
  # Extract data at points
  s <- as.data.frame(cbind(x = dp$x,y = dp$y, value = raster::extract(rdata, cbind(dp$x, dp$y))))
  depoint = s %>% st_as_sf(coords = c("x", "y"), dim = "XY") %>%
    st_set_crs(crsproj) %>% st_cast("MULTIPOINT")
  
  return(depoint)
}

# Generate areas
areasample <- function(n, rdata, crsproj){
  dearea <- aggregate(rdata, 100/n)
  #plot(dearea)
  # convert raster 2 to sp to be able to use extract()
  dearea <- rasterToPolygons(dearea)
  # Extract data at areas
  # sum values raster 1 in raster 2 with weights proportional to the overlapping area
  re <- raster::extract(rdata, dearea, weights = TRUE, normalizeWeights = TRUE)
  # this returns values and weights.
  # I cannot use extract() with fun = mean when weights = TRUE so I do it with values*weights
  dearea$value <- sapply(re, function(m){sum(m[, 1]*m[, 2], na.rm = TRUE)})
  
  dearea = st_as_sf(dearea)[,c(2,3)]
  st_crs(dearea) = crsproj
  
  return(dearea)
}


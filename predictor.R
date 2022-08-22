dearea  <- st_transform(dearea, 4326)
depoint <- st_transform(depoint, 4326)
library(raster)
library(sf)
library(rgeoboundaries)
library(ggplot2)
library(viridis)
library(stars)
library(ncdf4)
library(mapview)

dearea  <- st_transform(dearea, 4326)
depoint <- st_transform(depoint, 4326)

fnCreateRaster <- function(filename,savename,var){
  nc_data <- nc_open(filename)
  # Save the print(nc) dump to a text file
  {
    sink(savename)
    print(nc_data)
    sink()
  }
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- ncvar_get(nc_data, "time")
  value <- ncvar_get(nc_data, var) # store the data in a 3-dimensional array
  
  fillvalue <- ncatt_get(nc_data, var, "_FillValue")
  nc_close(nc_data)
  value[value == fillvalue$value] <- NA
  
  sum =  value[ , ,1]
  for(i in 2:12){
    sum <- sum + value[ , ,i]
  }
  annual_mean <- sum/12
  r <- raster(t(annual_mean), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
              crs=CRS("EPSG:4326"))
  return(r)
}

# temp
df <- read_ncdf("temp.nc" )
r_temp_2016 <- fnCreateRaster('temp.nc','tem_meata.txt','t2m')
area_bound <- st_as_sf(st_union(dearea))

temp_uk <- mask(r_temp_2016,area_bound)
plot(temp_uk)
spol <- rasterToPolygons(temp_uk, dissolve = F)
temp_sf <- st_as_sf(spol)
colnames(temp_sf) <- c('temp','geometry')

df <- read_ncdf("precipation.nc" )
r_precipation_2016 <- fnCreateRaster('precipation.nc','pre_meata.txt','tp')
precipation_uk <- mask(r_precipation_2016,area_bound)
plot(precipation_uk)
spol <- rasterToPolygons(precipation_uk, dissolve = F)
precipation_sf <- st_as_sf(spol)
colnames(precipation_sf) <- c('precipation','geometry')

de_new <- dearea
de_new$area <- as.numeric(st_area(de_new)/(1000^2))

de_new$group <- 1:169

dearea_combine <- de_new %>%
  st_join(precipation_sf)%>%
  group_by(group,value)%>%
  summarise(area_sum = sum(area),
            precipation = weighted.mean(precipation,area/sum(area)))

dearea_cov<- de_new%>%
  st_join(temp_sf)%>%
  group_by(group,value)%>%
  summarise(area_sum = sum(area),
            temp_mean = weighted.mean(temp,area/sum(area)))

dearea_combine$temp_mean <- dearea_cov$temp_mean-273.15

dearea_cov <- dearea_combine[,c(2,4,5,6)]

# point

depoint_cov <- depoint_cov[,c(1,6)] %>%
  st_join(temp_sf) %>%
  st_join(precipation_sf)
# road length 
# area
library(rlang)
library(osmdata)

fnHighwayL <- function(box){
  highways <- opq(bbox = box) %>%
    add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% osmdata_sf()
  if(!is.null(highways$osm_lines)){
    road <- sum(st_length(highways$osm_lines))
  } else{
    road <- 0
  }
  
  return(road)
}

road_length <- NULL
for (i in length(road_length)+1 : nrow(dearea_combine)){
  print(i)
  box <- unname(st_bbox(dearea_combine[i,]))
  road <- fnHighwayL(box)
  road_length <- c(road_length, road)
}

dearea_cov$road_length_1500 <- road_length

# buffered 
# Here use UTM and meters
buffer_margin_1500 = units::set_units(1500, m)
buffered_point_utm <- st_buffer(depoint, buffer_margin_1500)
buffered_point <- st_transform(buffered_point_utm,4326)
road_length <- NULL
for (i in length(road_length)+1 : nrow(buffered_point)){
  print(i)
  box <- unname(st_bbox(buffered_point[i,]))
  road <- fnHighwayL(box)
  road_length <- c(road_length, road)
}

depoint_new$road_length_1500 <- road_length

# Areal road proportion
#----
fnCreateAreaProp <- function(radius, route, area){
  
  utm_route = st_transform(route, 2158)
  buffer_margin = units::set_units(radius, m)
  road_prop <- utm_route%>%
    st_buffer(buffer_margin) %>% 
    st_transform(4326)%>%
    st_union() %>%
    st_area()/st_area(area)
  
  return(road_prop)
}
road_prop <- NULL

for (i in length(road_prop)+1 : nrow(dearea_combine)){
  print(i)
  box <- unname(st_bbox(dearea_combine[i,]))
  
  highways <- opq(bbox = box) %>%
    add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% 
    osmdata_sf()
  
  route <- highways$osm_lines
  if(!is.null(route)){prop <- fnCreateAreaProp(1000,route,dearea_combine[i,])} else {prop <- 0}
  
  road_prop<- c(road_prop, prop)
}
# point within buffered route
#----
fnFindRoadIn <-function(radius,box){
  highways <- opq(bbox = box) %>%
    add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% 
    osmdata_sf()
  buffer_margin = units::set_units(radius, m)
  if(!is.null(highways$osm_lines)){
    road_buffered <- highways$osm_lines%>%
      st_transform(4326) %>%
      st_buffer(buffer_margin) %>% 
      st_transform(4326)
    
    road_in_point <- as.numeric(nrow(st_intersection(depoint, road_buffered)))
  }else {road_in_point <- 0}
  
  return(road_in_point)
}
road_in <- NULL
point_in_area <- st_within(depoint,dearea)
for(i in length(road_in)+1 : nrow(point_in_area)){
  print(i)
  box <- unname(st_bbox(point_in_area[i,]))
  road_in_point <- fnFindRoadIn(radius,box)
  road_in <- c(road_in,road_in_point)
}

#----
st_write(depoint_cov, "depoint_cov.shp")
st_write(dearea_cov, "dearea_cov.shp")

de1 <-  st_read("depoint_cov.shp")%>%
  st_transform(crsproj)

de2 <- st_read("dearea_cov.shp")%>%
  st_transform(crsproj)

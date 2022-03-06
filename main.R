packages = c('raster','rgdal',
             'rnaturalearth',
             'viridis','rnaturalearthhires','sf','INLA','rgeoboundaries','tidyverse')

package.check = lapply(packages, FUN =function(x){
  if(!require(x ,character.only = T))
    install.packages(x)
  if(!(x %in% ( .packages()  ) ) )
    library(x ,character.only = T)})

source('mesh.R')
source('dwscaler.R')
source('tarConstr.R')
source('geosp.R')
source('meld.R')
bd.sf <- geoboundaries("United Kingdom")
pre.sf = bd.sf
cutoff = 0.25
max.edge = c(0.7, 0.7)
offset = c(-0.05, -0.05)

prior.range = c(2, 0.01)
prior.sigma = c(10, 0.01)

proN = 27700
pre.sf = bd.sf
cell.x = 50 
cell.y =50

# point data, load the data when wd is on data
p.df <- read.csv("2016PM2.5_avg.csv")

p.lon =p.df$longitude
p.lat = p.df$latitude
p.value = p.df$mean


# area data, load the data when wd is on data
str_name<-'gwr_pm25_2016.tif' 
glo_pm = raster(str_name)



rr <- mask(crop(glo_pm, bd.sf), bd.sf)
fa = 30
ra <- raster:: aggregate(rr, fact = fa, fun = mean)
spol<-rasterToPolygons(ra, dissolve = F)
area.sf =st_as_sf(spol)

mesh = Mesh(bd.sf,max.edge,cut.off = 0.25,offset)

plot(mesh)

area.pre = area.pre(pre.sf,cell.x,cell.y,proN =27700)

# the output contain two elements, one is estimated mean and confidence interval by col
# another is the same results but organized by row

a = dwscaler.spde (mesh,area.pre,
                         prior.range = c(2, 0.01),prior.sigma = c(10, 0.01),
                         p.lon,p.lat,p.value,area.sf,proN = 27700)

 
form = y ~ 0 + b0 + f(s, model = spde) 

b=geosp(form,mesh,p.lon,p.lat,p.value,area.pre,prior.range = c(2, 0.01),prior.sigma = c(10, 0.01))

c = meld(mesh,area.pre,
         p.lon,p.lat,p.value,area.sf,
         prior.range = c(2, 0.01),prior.sigma = c(10, 0.01),proN = 27700)

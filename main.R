packages = c('raster','rgdal',
             'rnaturalearth',
             'viridis','rnaturalearthhires','sf','INLA','rgeoboundaries','tidyverse')

package.check = lapply(packages, FUN =function(x){
  if(!require(x ,character.only = T))
    install.packages(x)
  if(!(x %in% ( .packages()  ) ) )
    library(x ,character.only = T)})

setwd("~/Documents/Project 1/resources/sptialM")
source('mesh.R')
source('check.r')
source('dwscalerv2.R')
source('tarConstr.R')
source('meld_v3.R')
source('prediction.R')

cutoff = 0.25
max.edge = c(0.7, 0.7)
offset = c(-0.05, -0.05)

prior.range = c(2, 0.01)
prior.sigma = c(10, 0.01)

proN = 4326
cell.x = 50 
cell.y =50

# point data, load the data when wd is on data
setwd("~/Documents/Project 1/data")
bd.sf = geoboundaries("United Kingdom")
st_transform(bd.sf, proN)
p.df = read.csv("2016PM2.5_avg.csv")
p.df = p.df[,c(2,3,4)]
colnames(p.df)<-c('value.p','py','px')

p.sf = p.df %>%
  st_as_sf(coords = c("px", "py"), dim = "XY") %>% 
  st_set_crs(proN) %>%
  st_cast("MULTIPOINT")

# area data, load the data when wd is on data

str_name<-'area_2016_pm2.5/gwr_pm25_2016.tif' 
glo_pm = raster(str_name)
rr <- mask(crop(glo_pm, bd.sf), bd.sf)
fa = 30
ra <- raster:: aggregate(rr, fact = fa, fun = mean)
spol<-rasterToPolygons(ra, dissolve = F)
area.sf = st_as_sf(spol)
st_transform(area.sf, proN)
colnames(area.sf) = c('avalue','geometry')
st_write(area.sf, paste0(getwd(), "/", "dearea.shp"), delete_layer = TRUE) 
ggplot(data = p.sf$geometry) + geom_sf()

mesh = Mesh(area.sf,offset =c(-0.05,-0.05), cut.off = 0.1)

plot(mesh)

#  require a continuous surface 

area.pre = target(pre.sf = area.sf, cell.x =cell.x,cell.y =  cell.y, proN = proN)

# if only want to predict point data, no need for target function 


a = dwscaler.spde (mesh = mesh, target = area.pre, areaPre.sf = NULL,
                   prior.range = c(2, 0.01), prior.sigma = c(10, 0.01),
                   p.sf = p.sf, area.sf = area.sf)

 
formula = y ~ 0 + b0 + f(s, model = spde) 

b_1 = meld(formula =formula, mesh = mesh,
         p.sf = p.sf, target = area.pre, 
         prior.range = c(2, 0.01), prior.sigma = c(10, 0.01))

b_2 = meld(formula =formula, 
           mesh = mesh, p.sf = NULL,
           area.sf = area.sf,
           target = area.pre,
           prior.range = c(2, 0.01),prior.sigma = c(10, 0.01))

c = meld(formula =formula,
         mesh = mesh, target = area.pre,
         p.sf = p.sf,area.sf = area.sf,
         prior.range = c(2, 0.01),prior.sigma = c(10, 0.01), proN = proN)

######### aggregation ##############

 fa = 100
ra <- raster:: aggregate(rr, fact = fa, fun = mean)
spol<-rasterToPolygons(ra, dissolve = F)
block = st_as_sf(spol)
st_transform(block, proN)
colnames(block) = c('value.a','geometry')

a_blcok = dwscaler.spde (mesh = mesh, target = area.pre, areaPre.sf = block,
                   prior.range = c(2, 0.01), prior.sigma = c(10, 0.01),
                   p.sf = p.sf, area.sf = area.sf)

# compare surface and block
plot(block)
group = block$geometry
plot(a)
pre.group = aggregate(a, group, mean)

plot(pre.group)

m = st_join(a, block,left = F)
n =  st_join(block, a,left = F)

plot(group)
coo = as.matrix(st_coordinates(a[,1]))
points(coo, col = "red")


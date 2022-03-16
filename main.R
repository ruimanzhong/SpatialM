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

bd.sf <- geoboundaries("United Kingdom")

pre.sf = bd.sf
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
colnames(area.sf) = c('value.a','geometry')

ggplot(data = p.sf$geometry) + geom_sf()

mesh = Mesh(bd.sf, max.edge, cut.off = 0.1, offset)

plot(mesh)

area.pre = area.pre(pre.sf, cell.x, cell.y, proN = proN)

# the output contain two elements, one is estimated mean and confidence interval by col
# another is the same results but organized by row

a = dwscaler.spde (mesh, area.pre, 
                   prior.range = c(2, 0.01), 
                   prior.sigma = c(10, 0.01), 
                   p.sf, area.sf)

 
formula = y ~ 0 + b0 + f(s, model = spde) 

b_1 = meld(formula =formula, mesh = mesh,
         p.sf = p.sf, area.pre = area.pre, 
         prior.range = c(2, 0.01), prior.sigma = c(10, 0.01))

b_2 = meld(formula =formula, 
           mesh = mesh, p.sf = NULL,
           area.sf = area.sf,
           area.pre = area.pre,
           prior.range = c(2, 0.01),prior.sigma = c(10, 0.01))

c = meld(formula =formula,
         mesh = mesh, area.pre = area.pre,
         p.sf = p.sf,area.sf = area.sf,
         prior.range = c(2, 0.01),prior.sigma = c(10, 0.01),proN = proN)

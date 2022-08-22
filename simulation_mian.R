library(raster)
library(RandomFields)
library(spatstat)
library(stars)
library(sf)
library(ggplot2)
library(INLA)
library(maptools)
library(ggpubr)

source('samplegenerator.r')
source("fnCreateMesh.R")
source('fnCheckInputsMelding.R')
source("fnCheckInputsDown.R")
source('fnPredictMelding.R')
source('fnPredictDown.R')

crsproj <- 4326

# Boundary region. Unit square
R <- ppp(x = runif(100), y = runif(100), c(0, 1), c(0, 1))
boundaryregion <- sf::st_as_sf(Window(R)) %>% sf::st_set_crs(crsproj)
plot(boundaryregion)

# Grid points for prediction. dpcontsurface
bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
x <- seq(bb[1] - 0.1, bb[3] + 0.1, length.out = 50)
y <- seq(bb[2] - 0.1, bb[4] + 0.1, length.out = 50)
coop <- expand.grid(x, y)
coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = crsproj)
dpcontsurface <- coop_sf %>% st_join(boundaryregion, left = FALSE)
ggplot() + geom_sf(data = dpcontsurface)
dppoint <- dpcontsurface

# Region study unit square
r <- raster(nrows = 100, ncols = 100, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
x <- coordinates(r)[, 1]
y <- coordinates(r)[, 2]

fnSimulate <- function(scenario, N, var, scale, nu, intercept, pnum, anum, unif){
  
  for(i in 1:N){
    print(i)
    
    # Simulate continuous surface
    model <- RMwhittle(scale = scale, nu = nu, var = var) + RMtrend(mean = intercept)
    simu <- RFsimulate(model = model, x = x, y = y)
    values(r) <- simu$variable1
    truesurface <- raster::extract(r, as.matrix(st_coordinates(dppoint)))
    
    # Sample observations in points and areas
    if(unif){p1 <- punifsample(pnum, r, crsproj)}
    if(!unif){p1 <- pimsample(pnum, r, crsproj)}
    a1 <- areasample(anum, r, crsproj)
    
    # Fit model with melding and downscaler approaches
    mesh <- fnCreateMesh(p1, boundaryregion)
    
    meld <- fnPredictMelding(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                             boundaryregion = boundaryregion, mesh = mesh)
    
    down <- fnPredictDown(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                          boundaryregion = boundaryregion, mesh = mesh)
    
    # Save truesurface
    write.csv(cbind(truesurface,meld[[1]]$pred_mean,down[[1]]$pred_mean), paste0("results/res/", scenario, "results", i, ".csv"), row.names = FALSE)
  
    
  }}

fnSimulateCov <- function(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif){
  
  for(i in 1:N){
    print(i)
    
    # Simulate continuous surface
    model <- RMwhittle(scale = scale,nu = nu, var = var ) + x * beta + RMtrend(mean = intercept)
    simu <- RFsimulate(model = model, x = x, y = y)
    values(r) <- simu$variable1
    truesurface <- raster::extract(r, as.matrix(st_coordinates(dppoint)))
    
    # Sample observations in points and areas
    if(unif){p1 <- punifsample(pnum, r, crsproj)}
    if(!unif){p1 <- pimsample(pnum, r, crsproj)}
    a1 <- areasample(anum, r, crsproj)
    
    # Fit model with melding and downscaler approaches
    mesh <- fnCreateMesh(p1, boundaryregion)
    
    meld <- fnPredictMeldingCov(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                             boundaryregion = boundaryregion, mesh = mesh)
    
    down <- fnPredictDownCov(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                          boundaryregion = boundaryregion, mesh = mesh)
    
    # Save truesurface
    write.csv(cbind(truesurface,meld[[1]]$pred_mean,down[[1]]$pred_mean), paste0("results/res/", scenario, "results", i, ".csv"), row.names = FALSE)
  }}

# Calculate MSE comparing true and estimated surfaces
fnCalculateMSE <- function(scenario, N){
  mse_meld <- NULL
  mse_down <- NULL
  for(i in 1:N){
    print(i)
    results <- read.csv(paste0("results/res/", scenario, "results", i, ".csv"))
    mse_meld <- c(mse_meld, Metrics::mse(results[ ,1], results[ ,2]))
    mse_down <- c(mse_down, Metrics::mse(results[ ,1], results[ ,3]))
  }
  write.csv(cbind(mse_meld, mse_down), paste0("results/mse/", scenario, "MSE", N, ".csv"), row.names = FALSE)
}


#####################################
# INI SIMULATIONS
#####################################

# Number iterations
N <- 100

# Common parameters
nu <- 2
intercept <- 10

# no_ocv
#----
#
# Scenario
sce <- "sce1"
var = 1; scale = 0.01; unif = T
for(pnum in c(10,50,100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce2"
var = 1; scale = 0.1; unif = T
for(pnum in c(10, 50, 100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce3"
var = 4; scale = 0.01; unif = T
for(pnum in c(10, 50, 100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce4"
var = 4; scale = 0.1; unif = T
for(pnum in c(10, 50, 100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce5"
var = 1; scale = 0.01; unif = F
for(pnum in c(10,50,100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce6"
var = 1; scale = 0.1; unif = F
for(pnum in c(50, 100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce7"
var = 4; scale = 0.01; unif = F
for(pnum in c(10, 50, 100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce8"
var = 4; scale = 0.1; unif = F
for(pnum in c(10, 50, 100)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulate(scenario, N, var, scale, nu, intercept, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# cov
#-----
# Scenario
sce <- "sce9"
var = 1; scale = 0.01; unif = T; beta = 2;
for(pnum in c(10,50,100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce10"
var = 1; scale = 0.1; unif = T;beta = 2;
for(pnum in c(10, 50, 100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

sce <- "sce11"
var = 4; scale = 0.01; unif = T;beta = 2;
for(pnum in c( 50, 100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce12"
var = 4; scale = 0.1; unif = T;beta = 2;
for(pnum in c(10, 50, 100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce13"
var = 1; scale = 0.01; unif = F;beta = 2;
for(pnum in c(10,50,100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce14"
var = 1; scale = 0.1; unif = F;beta = 2;
for(pnum in c(50)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    #fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce15"
var = 4; scale = 0.01; unif = F;beta = 2;
for(pnum in c(10, 50, 100)){
  for(anum in c(2)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce16"
var = 4; scale = 0.1; unif = F;beta = 2;
for(pnum in c(50)){
  for(anum in c(10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Results
#----
a = 2
final_1 = NULL
sce <- "sce13"
for(pnum in c(10, 50, 100)){
  for(anum in c(a)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    results <- read.csv(paste0("results/mse/", scenario, "MSE", N, ".csv"))
    
    final_1 <- rbind(final_1,
      data.frame(MSE = results[,1], Point = as.factor(pnum), Method = 'Melding'),
                   data.frame(MSE = results[,2], Point = as.factor(pnum), Method = 'Downscaler'))
    
  }}

final_2 = NULL
sce <- "sce14"
for(pnum in c(10, 50, 100)){
  for(anum in c(a)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    results <- read.csv(paste0("results/mse/", scenario, "MSE", N, ".csv"))
    
    final_2 <- rbind(final_2,
                     data.frame(MSE = results[,1], Point = as.factor(pnum), Method = 'Melding'),
                     data.frame(MSE = results[,2], Point = as.factor(pnum), Method = 'Downscaler'))
    
  }}

final_3 = NULL
sce <- "sce15"
for(pnum in c(10, 50, 100)){
  for(anum in c(a)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    results <- read.csv(paste0("results/mse/", scenario, "MSE", N, ".csv"))
    
    final_3 <- rbind(final_3,
                     data.frame(MSE = results[,1], Point = as.factor(pnum) , Method = 'Melding'),
                     data.frame(MSE = results[,2], Point = as.factor(pnum), Method = 'Downscaler'))
    
  }}

final_4 = NULL
sce <- "sce16"
for(pnum in c(10, 50, 100)){
  for(anum in c(a)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    results <- read.csv(paste0("results/mse/", scenario, "MSE", N, ".csv"))
    
    final_4 <- rbind(final_4,
                     data.frame(MSE = results[,1], Point = as.factor(pnum), Method = 'Melding'),
                     data.frame(MSE = results[,2], Point = as.factor(pnum), Method = 'Downscaler'))
    
  }}



p1 = ggplot(final_1, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot() +ggtitle('Scen5 var = 1 scale = 0.01 4 area')+ theme(text = element_text(size = 15))

p2 = ggplot(final_2, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot()+ggtitle('Scen6 var = 1 scale = 0.1 4 area')+ theme(text = element_text(size = 15))

p3 = ggplot(final_3, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot()+ggtitle('Scen7 var = 4 scale = 0.01 4 area')+ theme(text = element_text(size = 15))

p4 = ggplot(final_4, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot()+ggtitle('Scen8 var = 4 scale = 0.1 4 area')+ theme(text = element_text(size = 15))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE)

a= 10

p5 = ggplot(final_1, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot() +ggtitle('Scen5 var = 1 scale = 0.01 10 area')+ theme(text = element_text(size = 15))

p6 = ggplot(final_2, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot()+ggtitle('Scen6 var = 1 scale = 0.1 10 area')+ theme(text = element_text(size = 15))

p7 = ggplot(final_3, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot()+ggtitle('Scen7 var = 4 scale = 0.01 10 area')+ theme(text = element_text(size = 15))

p8 = ggplot(final_4, aes(x = Point, y= MSE, fill=Method)) +
  geom_boxplot()+ggtitle('Scen8 var = 4 scale = 0.1 10 area')+ theme(text = element_text(size = 15))

ggarrange(p5, p6, p7, p8, ncol = 2, nrow = 2, common.legend = TRUE)



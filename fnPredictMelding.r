#####################################################################
# INI DATA
#####################################################################

# 
# library(INLA)
# library(sp)
# library(rgdal)
# library(rgeoboundaries)
# library(geoR)
# 
# data(gambia)
# 
# # point data
# total <- aggregate(gambia$pos, by = list(gambia$x, gambia$y), FUN = length)
# d <- data.frame(x = total$Group.1,y = total$Group.2)
# sps <- SpatialPoints(d[, c("x", "y")], proj4string = CRS("+proj=utm +zone=28"))
# spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84"))
# d[, c("long", "lat")] <- coordinates(spst)
# depoint <- d[, c("long", "lat")]
# depoint$x <- depoint$long
# depoint$y <- depoint$lat
# depoint$value <- rnorm(nrow(depoint))
# 
# # areal data
# map <- geoboundaries("Gambia", "adm1")
# dearea <- as(map, "Spatial")
# dearea <- dearea[-c(7,8),] # polygons 7 and 8 give error
# dearea$value <- rnorm(nrow(dearea))

# # I put prediction points and areas the same as estimation
# dppoint <- depoint
# dparea <- dearea
# 
# # Plots. TODO: Functions to plot data and predictions
# plot(depoint$x, depoint$y)
# plot(dearea)
# 

#####################################################################
# END DATA
#####################################################################




#####################################################################
# INI MAIN FUNCTION
#####################################################################

# Predicts using the melding approach
# depoint dataset estimation point. Columns x, y, value
# dearea dataset estimation area. SpatialPolygonsDataFrame. polygon, value
# dppoint dataset prediction point
# dparea dataset prediction area
source('fncheck.R')
fnPredictMelding <- function(depoint, dearea, dppoint, dparea){

# Use 1 for points and 2 for areas
# datasets estimation
de1 <- depoint
de2 <- dearea
# datasets prediction
dp1 <- dppoint
dp2 <- dparea


# check I pass at least de1 or de2 to fit the model 
# check I pass at least dp1 or dp2 to predict
# check formats are OK
# check CRS points and areas are the same.

fnCheckInputsMelding(de1, de2, dp1, dp2)

# Logical values indicating what datasets I have
de1ToF <- !is.null(de1)
de2ToF <- !is.null(de2)
dp1ToF <- !is.null(dp1)
dp2ToF <- !is.null(dp2)

# Create mesh using estimation points. TODO: construct mesh better
locations <- cbind(de1$x, de1$y)
mesh <- fnCreateMesh(locations)

# Create spde and index
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)

# Projection matrices for points (estimation point and prediction point)
if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = cbind(de1$x, de1$y))}
if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = cbind(dp1$x, dp1$y))}

# Create projection matrix A for areas (estimation area and prediction area)
if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}

# Create stk.full
stk.full <- fnCreateStack(de1ToF, de2ToF, dp1ToF, dp2ToF)

# Specify formula melding
formula <- y ~ 0 + b0 + f(s, model = spde)

# Call inla()
res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))

# Retrieve predictions points
# Predictions points
if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, "pred1", dp1)}
# Predictions areas
if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, "pred2", dp2)}

}

#####################################################################
# END MAIN FUNCTION
#####################################################################



#####################################################################
# INI AUXILIARY FUNCTIONS
#####################################################################


# Check inputs are correct. TODO: this function
fnCheckInputsMelding <- function(de1, de2, dp1, dp2){
  
}

# Create mesh. TODO: default values for max.edge and cutoff
fnCreateMesh <- function(locations){
  mesh <- inla.mesh.2d(loc = locations, max.edge = c(0.1, 5), cutoff = 0.01)
  return(mesh)
}

# Create projection matrix A for areas. TODO: use sf instead of SpatialPolygonDataframe
fnProjectionMatrixArea <- function(dataset, mesh){
  spol <- dataset # SpatialPolygon
  spoints <- SpatialPoints(mesh$loc[, c(1,2)])
  crs(spoints) <- crs(spol)
  locin <- mesh$loc[, c(1,2)][as.vector(which(!is.na(over(spoints, spol)[,1]))),]
  # over(spoints, spol)[,1]. recupero primera columna
  block <- rep(0, nrow(locin))
  for(i in 1:length(spol)){
    spoints <- SpatialPoints(locin)
    crs(spoints) <- crs(spol)
    block[as.vector(which(!is.na(over(spoints, spol[i]))))] <- i
  }
  A <- inla.spde.make.A(mesh = mesh, loc = locin, block = block, block.rescale = "sum")
  return(A)
}

# Create stk.full
fnCreateStack <- function(de1ToF, de2ToF, dp1ToF, dp2ToF){
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value), A = list(1, Ae1), effects = list(data.frame(b0 = rep(1, nrow(de1))), s = indexs))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value), A = list(1, Ae2), effects = list(data.frame(b0 = rep(1, nrow(de2))), s = indexs))}
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA), A = list(1, Ap1), effects = list(data.frame(b0 = rep(1, nrow(dp1))), s = indexs))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(1, Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), s = indexs))}
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  return(stk.full)
}


# Retrieve predictions
fnRetrievePredictions <- function(stack, tag, dataset){
  index <- inla.stack.index(stack = stack, tag = tag)$data
  dataset$pred_mean <- res$summary.fitted.values[index, "mean"]
  dataset$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  dataset$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  return(dataset)
}

#####################################################################
# END AUXILIARY FUNCTIONS
#####################################################################

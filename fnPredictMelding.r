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
source('fnCreateMesh.R')
fnPredictMelding <- function(depoint = NULL, dearea = NULL, 
                             dppoint = NULL, dparea = NULL, bd, mesh = NULL,
                             prior.range = NULL, prior.sigma = NULL){

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

fnCheckInputsMelding(de1, de2, dp1, dp2, bd)

# Logical values indicating what datasets I have
de1ToF <- !is.null(de1)
de2ToF <- !is.null(de2)
dp1ToF <- !is.null(dp1)
dp2ToF <- !is.null(dp2)
meshToF <- is.null(mesh)
spdeIndex <- is.null(prior.range)
# Create mesh using estimation points. TODO: construct mesh better
#locations <- cbind(de1$x, de1$y)

if(meshToF){mesh <- fnCreateMesh(de1,de2,bd)}


# Create spde and index
if(spdeIndex){
  spde = inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
} else {
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = prior.range, 
    prior.sigma = prior.sigma)
}

indexs <- inla.spde.make.index("s", spde$n.spde)

# Projection matrices for points (estimation point and prediction point)
if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}

# Create projection matrix A for areas (estimation area and prediction area)
if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}

# Create stk.full
message('Creating stack')
stk.e1 = NULL
stk.e2 = NULL
stk.p1 = NULL
stk.p2 = NULL

# estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$pvalue), A = list(1, Ae1), effects = list(data.frame(b0 = rep(1, nrow(de1))), s = indexs$s))}
if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$avalue), A = list(1, Ae2), effects = list(data.frame(b0 = rep(1, nrow(de2))), s = indexs$s))}
if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA), A = list(1, Ap1), effects = list(data.frame(b0 = rep(1, nrow(dp1))), s = indexs$s))}
if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(1, Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), s = indexs$s))}
# construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])


# Specify formula melding
formula <- y ~ 0 + b0 + f(s, model = spde)

# Call inla()
res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))

# Retrieve predictions points
# Predictions points
if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, res, "pred1", dp1)}
# Predictions areas
if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, res, "pred2", dp2)}

return(list(dp1,dp2))
}

#####################################################################
# END MAIN FUNCTION
#####################################################################



#####################################################################
# INI AUXILIARY FUNCTIONS
#####################################################################


# Check inputs are correct. TODO: this function
# fnCheckInputsMelding <- function(de1, de2, dp1, dp2){
#   
# }

# Create mesh. TODO: default values for max.edge and cutoff
# fnCreateMesh <- function(locations){
#   mesh <- inla.mesh.2d(loc = locations, max.edge = c(0.1, 5), cutoff = 0.01)
#   return(mesh)
# }

# Create projection matrix A for areas.
fnProjectionMatrixArea <- function(de2, mesh){
  message('Creating areal projection matrix')
  meshcoo = data.frame(X = mesh$loc[,1], Y = mesh$loc[,2])
  
  meshin = meshcoo %>%
    st_as_sf(coords = c("X", "Y"), dim = "XY") %>% 
    st_set_crs(proN) %>%
    st_cast("MULTIPOINT")
  
  # find points in mesh n area.sf
  locin= st_join(meshin, de2, left = F)
  
  block <- rep(0, nrow(locin))
  for(i in 1:nrow(de2)) {
    block[as.vector(which(!is.na(st_join(locin, de2[i,], left = T)$avalue.y)))] <- i
  } 
  
  A <- inla.spde.make.A(mesh=mesh, 
                        loc=as.matrix(st_coordinates(locin[,1]))[,c(1,2)], 
                        block=block, 
                        block.rescale="sum")
  return(A)
}

# Create stk.full
# fnCreateStack <- function(de1ToF, de2ToF, dp1ToF, dp2ToF, de1, de2, dp1, dp2){
#   message('Creating stack')
#    stk.e1 = NULL
#    stk.e2 = NULL
#    stk.p1 = NULL
#    stk.p2 = NULL
#    
#    # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
#   if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$pvalue), A = list(1, Ae1), effects = list(data.frame(b0 = rep(1, nrow(de1))), s = indexs$s))}
#   if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$avalue), A = list(1, Ae2), effects = list(data.frame(b0 = rep(1, nrow(de2))), s = indexs$s))}
#   if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA), A = list(1, Ap1), effects = list(data.frame(b0 = rep(1, nrow(dp1))), s = indexs$s))}
#   if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(1, Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), s = indexs$s))}
#   # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
#   stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
#   return(stk.full)
# }


# Retrieve predictions
fnRetrievePredictions <- function(stack, res, tag, dataset){
  index <- inla.stack.index(stack = stack, tag = tag)$data
  dataset$pred_mean <- res$summary.fitted.values[index, "mean"]
  dataset$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  dataset$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  return(dataset)
}

#####################################################################
# END AUXILIARY FUNCTIONS
#####################################################################

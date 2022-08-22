#' @param prior.sigma A two-dimensional vector to construct SPDE. More information can be found inla.spde2.pcmatern()
#' @param depoint A sf obj includes 'pvalue' and 'geometry' columns, the point data
#' @param dearea A sf obj includes 'avalue' and 'geometry' columns, the areal data
#' @return A dataframe with predicted mean, 95% credible interval of the prediction

fnPredictDown <- function(depoint, dearea, dppoint = NULL, dparea = NULL, boundaryregion,
                          mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  
  
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Set avalue and pvalue instead of value
  de1$pvalue <- de1$value
  de2$avalue <- de2$value
  dp1$pvalue <- dp1$value
  dp2$avalue <- dp2$value  
  
  # Check inputs
  
  # check function
  fnCheckInputsDown(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  
  # Mesh
  if (is.null(mesh) == T) {
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  indexs <- inla.spde.make.index("s", spde$n.spde)
  indexs1 <- inla.spde.make.index("s1", spde$n.spde)
  
  
  
  # Match Areal data and point data, Match dppoint and dearea
  locin <- st_join(de1, de2, left = F)
  
  # If TRUE, predict in points
  if(dp1ToF){
    locin_pred <- st_join(dp1, de2, left = FALSE)
  }
  
  # If TRUE, predict in areas. Construct prediction surface
  if(dp2ToF){
    # Bounding box
    bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
    # Grid
    x <- seq(bb[1] - 1, bb[3] + 1, length.out = 50)
    y <- seq(bb[2] - 1, bb[4] + 1, length.out = 50)
    coop <- expand.grid(x, y)
    coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = st_crs(de1))
    # Keep points inside
    dpcontsurface <- coop_sf %>% st_join(boundaryregion, left = FALSE) %>% st_join(de2, left = FALSE)
  }
  
  
  
  # Create projection matrices A
  # Matrix A estimation
  Ae <- inla.spde.make.A(mesh, loc = as.matrix(st_coordinates(locin[, 1]))[, c(1, 2)])
  # Matrix A points prediction
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(locin_pred)[, c(1,2)]))}
  # Matrix A areas prediction (need to predict in a continuous surface first, dense points)
  if(dp2ToF){Ap2 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dpcontsurface[, 1])))}
  
  # stack estimation
  stk.p1 <- NULL
  stk.p2 <- NULL
  
  stk.e <- inla.stack(tag = "est", data = list(y = locin$pvalue),
                      A = list(1, Ae * locin$avalue, Ae),
                      effects = list(data.frame(b0 = rep(1, nrow(locin)), X = locin$avalue), s = indexs$s, s1 = indexs1$s1))
  # stack prediction points
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),
                                  A = list(1, Ap1 * locin_pred$avalue, Ap1),
                                  effects = list(data.frame(b0 = rep(1, nrow(locin_pred)), X = locin_pred$avalue), s = indexs$s, s1 = indexs1$s1))}
  # stack prediction areas
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA),
                                  A = list(1, Ap2 * dpcontsurface$avalue, Ap2),
                                  effects = list(data.frame(b0 = rep(1, nrow(dpcontsurface)), X = dpcontsurface$avalue), s = indexs$s, s1 = indexs1$s1))}
  # construct stack full with the data we have
  stk.full <-  do.call(inla.stack, list(stk.e, stk.p1, stk.p2)[c( T , dp1ToF, dp2ToF)])
  
  
  # Specify formula downscaling
  formula <- y ~ 0 + b0 + X + f(s, model = spde) + f(s1, model = spde)
  
  
  # Call inla()
  res <- inla(formula, data = inla.stack.data(stk.full), control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)))
  
  
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions_p(stk.full, res, "pred1", locin_pred)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions_a(stk.full, res, "pred2", dpcontsurface, dp2)}
  
  return(list(dp1, dp2))
}


#####################################################################
# INI AUXILIARY FUNCTIONS
#####################################################################

# Retrieve predictions
fnRetrievePredictions_p <- function(stack, res, tag, dataset){
  index <- inla.stack.index(stack = stack, tag = tag)$data
  dataset$pred_mean <- res$summary.fitted.values[index, "mean"]
  dataset$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  dataset$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  
  return(dataset)
}

fnRetrievePredictions_a <- function(stack, res, tag, dataset, dp){
  index <- inla.stack.index(stack = stack, tag = tag)$data
  dataset$pred_mean <- res$summary.fitted.values[index, "mean"]
  dataset$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  dataset$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  agg <- aggregate(dataset, dp[2], mean) # PAULA: This should come from inla.posterior.sample
  return(agg)
}

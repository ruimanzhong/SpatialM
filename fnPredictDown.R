#' @param prior.sigma A two-dimensional vector to construct SPDE. More information can be found inla.spde2.pcmatern()
#' @param depoint A sf obj includes 'pvalue' and 'geometry' columns, the point data
#' @param dearea A sf obj includes 'avalue' and 'geometry' columns, the areal data
#' @return A dataframe with predicted mean, 95% credible interval of the prediction

source("fnCreateMesh.R")
source("fnCheckInputsDown.R")
fnPredictDown <- function(depoint, dearea, dppoint = NULL, dparea = NULL, boundaryregion
                          , mesh = NULL, priorspdesigma = NULL, priorspderange = NULL) {
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # check function
  # TO DO : check crs of input is crsproj
  fnCheckInputsDown(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  # Match Areal data and point data, Match dppoint and dearea
  locin <- st_join(de1, de2, left = F)
  
  # Projection Matrix for point data and dppoint area
  
  if (is.null(mesh) == T) {
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else {
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  
  # set formula
  formula <- y ~ 0 + b0 + X + f(s, model = spde) + f(s1, model = spde)
  indexs <- inla.spde.make.index("s", spde$n.spde)
  indexs1 <- inla.spde.make.index("s1", spde$n.spde)
  
  if(dp1ToF){
  locin_pred <- st_join(dp1,de2, left = FALSE)
  }
  
  if(dp2ToF){
    # Construct Prediction Surface
    bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
    
    # Grid
    x <- seq(bb[1] - 1, bb[3] + 1, length.out = 50)
    y <- seq(bb[2] - 1, bb[4] + 1, length.out = 50)
    coop <- expand.grid(x, y)
    coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = st_crs(de1))
    dpcontsurface <- coop_sf %>% st_join(boundaryregion, left = FALSE) %>% st_join(de2, left = FALSE)
  }

  # construct stack
  Ae <- inla.spde.make.A(mesh, loc = as.matrix(st_coordinates(locin[, 1]))[, c(1, 2)])
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(locin_pred )[ , c(1,2)]))}
  if(dp2ToF){Ap2 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dpcontsurface[, 1])))}
  
  
  stk.e <- inla.stack(tag = "est", data = list(y = locin$pvalue), A = list(1, Ae * locin$avalue, Ae), effects = list(data.frame(b0 = rep(1, nrow(locin)), X = locin$avalue), s = indexs$s, s1 = indexs1$s1))
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA), A = list(1, Ap1 * locin_pred$avalue, Ap1), effects = list(data.frame(b0 = rep(1, nrow(locin_pred)), X = locin_pred$avalue), s = indexs$s, s1 = indexs1$s1))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(1, Ap2 * dpcontsurface$avalue, Ap2), effects = list(data.frame(b0 = rep(1, nrow(dpcontsurface)), X = dpcontsurface$avalue), s = indexs$s, s1 = indexs1$s1))}
  
  stk.full <-  do.call(inla.stack, list(stk.e, stk.p1, stk.p2)[c( T , dp1ToF, dp2ToF)])
  
  # estimation and prediction
  res <- inla(formula,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              )
  )
  
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions_p(stk.full, res, "pred1", locin_pred)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions_a(stk.full, res, "pred2",dpcontsurface, dp2)}
  
  return(list(dp1,dp2))
}

# ----
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
  res = aggregate(dataset, dp[2], mean)
  return(res)
}
# Prior Check
fnCheckPrior <- function(prior.range, prior.sigma){
  if(length(prior.range) != 2 ||sum(prior.range >0) != 2 || prior.range[2] > 1){
    stop("The parameters of prior.sigma are not valid, the length 2 vector contains spatial range (>0) and 
         prior probability of the deviation (0~1)")
  }
  if(length(prior.sigma) != 2 || sum(prior.sigma >0) != 2 || prior.sigma[2] > 1){
    stop("The parameters of prior.sigma are not valid, the length 2 vector contains marginal standard deviation (>0) and 
         prior probability of the deviation (0~1)")
  }
}

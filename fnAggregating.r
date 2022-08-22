fnInla <- function(depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                          mesh = NULL, priorspdesigma = NULL, priorspderange = NULL,N){
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Check inputs
  fnCheckInputsMelding(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  de1ToF <- !is.null(de1)
  de2ToF <- !is.null(de2)
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  
  # Mesh
  if(is.null(mesh)){
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
  
  
  # Projection matrices for points (estimation point and prediction point)
  if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}
  
  # Create projection matrix A for areas (estimation area and prediction area)
  if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
  if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}
  
  # Create stk.full
  stk.e1 = NULL
  stk.e2 = NULL
  stk.p1 = NULL
  stk.p2 = NULL
  
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value), A = list(1, Ae1), effects = list(data.frame(b0 = rep(1, nrow(de1))), s = indexs$s))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value), A = list(1, Ae2), effects = list(data.frame(b0 = rep(1, nrow(de2))), s = indexs$s))}
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA), A = list(1, Ap1), effects = list(data.frame(b0 = rep(1, nrow(dp1))), s = indexs$s))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(1, Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), s = indexs$s))}
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  
  
  
  # Specify formula melding
  formula <- y ~ 0 + b0 + f(s, model = spde)
  
  # Call inla()
  res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
              control.compute=list(dic=TRUE, config = TRUE),
              control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
  samples <- inla.posterior.sample(N,res)
  predicted_samples<- lapply(samples,function(x) x$latent[as.numeric(nrow(stk.e1$A)+nrow(stk.e2$A)+1):nrow(stk.full$A)])
  
  return(predicted_samples)
}

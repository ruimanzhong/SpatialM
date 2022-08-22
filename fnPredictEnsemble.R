library(caret)
library(glmnet)
source("fnPredictDown.R")
source("fnPredictMelding.R")
#source('CV.R')
fnPredictEnsemble = function(depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                             mesh = NULL, priorspdesigma = NULL, priorspderange = NULL, SCV = F ){
  
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  prior1 <- priorspdesigma
  prior2 <- priorspderange
  mesh <- mesh
  set.seed(123)
  # Logical values indicating what datasets I have
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  dt = sort(sample(nrow(de1), nrow(de1)*.8))
  
  
  de_train <- de1[dt,]
  de_test <- de1[-dt,]
  
  base_melding <- fnPredictMelding(depoint = de_train, de2, dppoint = de_test, dparea = NULL, boundaryregion,
                                   mesh, priorspdesigma = prior1, priorspderange = prior2)[[1]]
  
  base_down <- fnPredictDown(depoint = de_train, de2, dppoint = de_test, dparea = NULL, boundaryregion,
                             mesh, priorspdesigma = prior1, priorspderange = prior2)[[1]]
  data <- st_join(base_melding,base_down)
  
  final_melding <- fnPredictMelding(depoint = de1, de2, dppoint = dp1, dparea = dp2, boundaryregion,
                                    mesh, priorspdesigma = prior1, priorspderange = prior2)[[1]]
  
  final_down <- fnPredictDown(depoint = de1, de2, dppoint = dp1, dparea = dp2, boundaryregion,
                              mesh, priorspdesigma = prior1, priorspderange = prior2)[[1]]
  
  newdata <- st_join(final_melding,final_down)
  
  #standard CV 
  #----
  if(SCV == F) {
    regressControl  <- trainControl(method="repeatedcv",
                                    number = 3,
                                    repeats = 3
    ) 
    stack_regression_const <- train(
      value ~ pred_mean.x + pred_mean.y, data = data, method = "lm",
      tuneGrid  = expand.grid(intercept = FALSE),
      trControl = regressControl,
      tuneLength = 10
    )
    
    final_pred = predict(stack_regression_const, newdata = newdata)
  }
  #----
  # spatial coefficients
  #----
  # Specify formula 
  if(SCV == T){
    
    mesh_ensem <- fnCreateMesh(data, boundaryregion)
    
    formula <- y ~ 0  + X1 + X2 + f(s, model = spde) + f(s1, model = spde)
    # 
    if(!is.null(prior1) & !is.null(prior2)){
      fnCheckPrior(prior1, prior2)
      spde <- inla.spde2.pcmatern(mesh = mesh_ensem, prior.range = prior2, prior.sigma = prior1)
    }else{
      spde <- inla.spde2.matern(mesh = mesh_ensem, alpha = 2, constr = T)
    }
    
    indexs <- inla.spde.make.index("s", spde$n.spde)
    indexs1 <- inla.spde.make.index("s1", spde$n.spde)
    
    Ae <- inla.spde.make.A(mesh_ensem, loc = as.matrix(st_coordinates(data$geometry))[, c(1, 2)])
    Ap <- inla.spde.make.A(mesh = mesh_ensem, loc = as.matrix(st_coordinates(dp1$geometry)))
    # 
    stk.e <- inla.stack(tag = "est", data = list(y = data$value),
                        A = list( Ae * data$pred_mean.x, Ae * data$pred_mean.y, Ae),
                        effects = list(data.frame(X1 = data$pred_mean.x, X2 = data$pred_mean.y), s = indexs$s, s1 = indexs1$s1))
    stk.p <- inla.stack(tag = "pred2", data = list(y = NA),
                        A = list( Ap * data$pred_mean.x, Ap * data$pred_mean.y, Ap ),
                        effects = list(data.frame( X1 = data$pred_mean.x, X2 = data$pred_mean.y), s = indexs$s, s1 = indexs1$s1))
    # # Call inla()
    res <- inla(formula, data = inla.stack.data(stk.full), control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)))
    
  } 
  
  
  
  return(final_pred)
}

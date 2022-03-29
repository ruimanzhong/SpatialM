geosp = function(formula, mesh, p.sf, target, areaPre.sf = NULL, spde, proN = 4326){
  # Projection Matrix
  Ap = inla.spde.make.A(mesh , as.matrix(st_coordinates(p.sf[,1]))[,c(1,2)])
  Apre = inla.spde.make.A(mesh = mesh, loc = as.matrix (st_coordinates(target[,1])))
  
  # construct stack
  indexs = inla.spde.make.index("s", spde$n.spde)
  
  #construct stack
  stk.full = stack.full.geo(Ap,Apre,p.sf,target,indexs)
  
  #estimation and prediction
  res <- inla(formula,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              ))
  
  #prediction data frame
return(pred(res,stk.full,target,areaPre.sf,proN)) 
}

#*************Construct stack*********************
#*
stack.full.geo = function(Ap,Apre,p.sf,target,indexs){
  stk.e <- inla.stack(
  tag = "est",
  data = list(y = p.sf$value),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(p.sf))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Apre),
  effects = list(data.frame(b0 = rep(1, nrow(st_coordinates(target[,1])))), s = indexs)
)

return (inla.stack(stk.e, stk.p))

}

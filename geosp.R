geosp = function(form,mesh,p.lon,p.lat,p.value,area.pre,
                 prior.range = c(2, 0.01),prior.sigma = c(10, 0.01)){
  
  p.df = data.frame(lon = p.lon,lat = p.lat,value = p.value )
  
  # Projection Matrix
  A = inla.spde.make.A(mesh, cbind(p.lon,p.lat))
  Ap = inla.spde.make.A(mesh = mesh, loc = as.matrix (st_coordinates(area.pre[,1])))
  
  #construct SPDE
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = prior.range, 
    prior.sigma = prior.sigma)
  
  # construct stack
  indexs = inla.spde.make.index("s", spde$n.spde)
  
  #construct stack
  stk.full = stack.full.geo(A,Ap,p.df,area.pre,indexs)
  
  #estimation and prediction
  res <- inla(form,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              ))
  
  #prediction data frame
 pred(res,stk.full,area.pre)
}

#*************Construct stack*********************
#*
stack.full.geo = function(A,Ap,p.df,area.pre,indexs){
stk.e <- inla.stack(
  tag = "est",
  data = list(y = p.df$value),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(p.df))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(st_coordinates(area.pre[,1])))), s = indexs)
)

inla.stack(stk.e, stk.p)

}

meld = function(mesh,area.pre,
                p.lon,p.lat,p.value,area.sf,
                prior.range = c(2, 0.01),prior.sigma = c(10, 0.01),proN = 27700){
  
  st_crs(area.sf) <- proN
  p.df = data.frame(lon = p.lon,lat = p.lat,value = p.value )
  colnames(area.sf) = c("value.a", "geometry")
  attributes(area.sf)
  
  locin_pred = st_join(area.pre,area.sf,left = F)
  
  # Projection Matrix
  meshcoo = data.frame(long = mesh$loc[,1],lat = mesh$loc[,2] )
  
  meshin = meshcoo %>%
    st_as_sf(coords = c("long", "lat"), dim = "XY") %>% 
    st_set_crs(27700) %>%
    st_cast("MULTIPOINT")
  
  # find points in mesh n area.sf
  locin= st_join(meshin,area.sf,left = F)
  # find which block each point lies
  
  block <- rep(0, nrow(locin))
  
  for(i in 1:nrow(area.sf)){
    block[as.vector(which(!is.na(st_join(locin, area.sf[i,],left = T)$value.a.y)))] <- i
  } 
  
  Ap = inla.spde.make.A(mesh, cbind(p.lon,p.lat))
  
  Aa <- inla.spde.make.A(mesh=mesh, 
                         loc=as.matrix(st_coordinates(locin[,1]))[,c(1,2)], 
                         block=block, 
                         block.rescale="sum")
  
  Apred = inla.spde.make.A(mesh = mesh, loc = as.matrix (st_coordinates(area.pre[,1])))
  
  #construct SPDE
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = prior.range, 
    prior.sigma = prior.sigma) 
  
  # construct stack

  stk.full = stack.full.meld(Ap,Aa,Apred,p.df,area.sf,locin_pred,spde)
  
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

#******************Construct Stack*********************************
stack.full.meld = function(Ap,Aa,Apred,p.df,area.sf,locin_pred,spde){
  ya = area.sf$value.a
  
  stk.p <- inla.stack(tag='point',
                      data=list(y=p.df$value),
                      A=list(Ap, 1),effects=list(s=1:spde$n.spde, data.frame(b0 =rep(1, length(p.df$value)))))
  
  stk.a <- inla.stack(tag= 'areal',data=list(y=ya),
                      A=list(Aa, 1),effects=list(s=1:spde$n.spde, 
                                                 data.frame(b0 =rep(1, length(ya)))))
  
  stk.pred <- inla.stack(tag= 'pred',
                         data=list(y=NA),
                         A=list(Apred, 1),
                         effects=list(s=1:spde$n.spde, 
                                      data.frame(b0 =rep(1, nrow(locin_pred)))))
  
 inla.stack(stk.p, stk.a, stk.pred)
}


dwscaler.spde = function(mesh,area.pre,
                         prior.range = c(2, 0.01),prior.sigma = c(10, 0.01),
                         p.lon,p.lat,p.value,area.sf,proN = 27700){
   #RENAME area.sf
  colnames(area.sf) = c("value.a", "geometry")
  st_crs(area.sf) <- proN
   # Match Areal data and point data, Match area.pre and area.sf
  locin = Match(p.lon,p.lat,p.value,area.sf)
  locin_pred = st_join(area.pre,area.sf,left = F)
  
  #Projection Matrix for point data and target area
  A = inla.spde.make.A(mesh, cbind(p.lon,p.lat))
  Ap = inla.spde.make.A(mesh = mesh, loc = as.matrix (st_coordinates(area.pre[,1])))
  
  #construct SPDE
  spde = inla.spde2.pcmatern(
     mesh = mesh,
     prior.range = prior.range, 
     prior.sigma = prior.sigma) 
  
  # set formula
  form <- y ~ 0 + b0  + X + f(s, model = spde)+f(s1,model= spde)
  indexs = inla.spde.make.index("s", spde$n.spde)
  indexs1 = inla.spde.make.index("s1", spde$n.spde)
  
  #construct stack
  stk.full = stack.full.ds(locin,A,Ap,locin_pred,indexs,indexs1)
  
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

#***************Match Areal data and point data ************

Match = function(p.lon,p.lat,p.value,area.sf,proN = 27700 ){
  area=area.sf %>% 
    st_set_crs(proN)
  
  p.df = data.frame(lon = p.lon,lat = p.lat,y = p.value )
  p.df %>%
    st_as_sf(coords = c("lon", "lat"), dim = "XY") %>% st_set_crs(proN) %>%
    st_cast("MULTIPOINT")%>%
    st_join(area,left = F)
}


#*************Construct stack*************
stack.full.ds = function(locin,A,Ap,locin_pred,indexs,indexs1){
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = locin$y),
    A = list(1,A*locin$value.a,A),
    effects = list(data.frame(b0 = rep(1, nrow(locin)),X = locin$value.a), s = indexs, s1 = indexs1)
  )
  
  # stack for prediction stk.p
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap*locin_pred$value.a,Ap),
    effects = list(data.frame(b0 = rep(1, nrow(locin_pred)), X= locin_pred$value.a), s = indexs, s1 = indexs1)
  )
  
  inla.stack(stk.e, stk.p)
}


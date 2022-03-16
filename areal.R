areal = function(formula, area.sf, area.pre, mesh, spde, proN) {

  
  locin_pred = st_join(area.pre, area.sf, left = F)
  
  # Projection Matrix
  meshcoo = data.frame(long = mesh$loc[,1], lat = mesh$loc[,2] )
  
  meshin = meshcoo %>%
    st_as_sf(coords = c("long", "lat"), dim = "XY") %>% 
    st_set_crs(proN) %>%
    st_cast("MULTIPOINT")
  
  # find points in mesh n area.sf 
  locin = st_join(meshin, area.sf, left = F)
  
  Aa = inla.spde.make.A (mesh = mesh, 
                         loc = as.matrix(st_coordinates(locin[,1]))[,c(1,2)]
                         )
  
  Apred = inla.spde.make.A (mesh = mesh, 
                            loc = as.matrix (st_coordinates(area.pre[,1]) 
                                             ) 
                            )
  
  stk.full = stack.full.areal(Aa,Apred,area.sf,locin_pred,spde)
  
  res = inla ( formula,
               data = inla.stack.data(stk.full),
               control.predictor = list( 
                 compute = TRUE, 
                 A = inla.stack.A(stk.full)
                 )
             )
  
  return(pred(res, stk.full, area.pre))
}

stack.full.areal = function(Aa, Apred, area.sf, locin_pred, spde) {
  
  stk.a = inla.stack( tag= 'areal',
                      data=list(y=locin$value.a),
                      A=list(Aa, 1),effects=list( s=1:spde$n.spde, 
                                                  data.frame(b0 =rep(1, length(locin$value.a)
                                                                    )
                                                            )
                                                 )
                      )
  
  stk.pred = inla.stack( tag= 'pred',
                          data=list(y=NA),
                          A=list(Apred, 1),
                          effects=list(s=1:spde$n.spde, 
                                       data.frame(b0 =rep(1, nrow(locin_pred)
                                                         )
                                                 )
                                      )
                         )
  
  return(inla.stack(stk.a, stk.pred))
}
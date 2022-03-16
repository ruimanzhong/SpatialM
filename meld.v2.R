source('check.r')
source('mesh.R')
meld = function(mesh = NULL,area.pre,formula,
                p.sf = NULL ,area.sf = NULL,
                prior.range = c(2, 0.01),prior.sigma = c(10, 0.01),proN = 27700){
  
  check_meld(area.sf, p.sf,area.pre,proN)
  
  Ap  =NULL
  Aa = NULL
  
  if(is.null(mesh) ==T){
    message("Generate mesh by default")
    
    mesh = Mesh(bd.sf,max.edge = c(0.7, 0.7),
                cut.off = 0.25,offset = c(-0.05, -0.05))
  }
  
  if(is.null(area.sf) == F){
    
    message("Join the point data and areal data")
    locin_pred = st_join(area.pre,area.sf,left = F)
    
    # Projection Matrix
    meshcoo = data.frame(long = mesh$loc[,1],lat = mesh$loc[,2] )
    
    meshin = meshcoo %>%
      st_as_sf(coords = c("long", "lat"), dim = "XY") %>% 
      st_set_crs(proN) %>%
      st_cast("MULTIPOINT")
    
    # find points in mesh n area.sf
    locin= st_join(meshin,area.sf,left = F)
    
    
    message('Constructing Projection Matrix')
    
    block <- rep(0, nrow(locin))
    
    for(i in 1:nrow(area.sf)){
      block[as.vector(which(!is.na(st_join(locin, area.sf[i,],left = T)$value.a.y)))] <- i
    } 
    
    Aa <- inla.spde.make.A(mesh=mesh, 
                           loc=as.matrix(st_coordinates(locin[,1]))[,c(1,2)], 
                           block=block, 
                           block.rescale="sum")
  }
  
  if(is.null(p.sf) == F){
    Ap = inla.spde.make.A(mesh, as.matrix(st_coordinates(p.sf[,1]))[,c(1,2)])
  }
 
  
  
  Apred = inla.spde.make.A(mesh = mesh, loc = as.matrix (st_coordinates(area.pre[,1])))
  
  #construct SPDE
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = prior.range, 
    prior.sigma = prior.sigma) 
  
  # construct stack
  check_projM(Aa,Ap)
  
  stk.full = stack.full.meld(Ap,Aa,Apred,p.df,area.sf,locin_pred,spde)
  
  #estimation and prediction
  res <- inla(formula,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              ))
  
  #prediction data frame
  return(pred(res,stk.full,area.pre))
}

#******************Construct Stack*********************************
stack.full.meld = function(Ap,Aa,Apred,p.sf,area.sf,locin_pred,spde){
 
  stk.p <- inla.stack(tag='point',
                      data=list(y=p.sf$value.p),
                      A=list(Ap, 1),effects=list(s=1:spde$n.spde, data.frame(b0 =rep(1, length(p.sf$value.p)))))
  
  stk.a <- inla.stack(tag= 'areal',data=list(y=area.sf$value.a),
                      A=list(Aa, 1),effects=list(s=1:spde$n.spde, 
                                                 data.frame(b0 =rep(1, length(ya)))))
  
  stk.pred <- inla.stack(tag= 'pred',
                         data=list(y=NA),
                         A=list(Apred, 1),
                         effects=list(s=1:spde$n.spde, 
                                      data.frame(b0 =rep(1, nrow(locin_pred)))))
  
  inla.stack(stk.p, stk.a, stk.pred)
}

stack.full.geo = function(Ap,Apre,p.sf,area.pre,indexs){
  
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
    effects = list(data.frame(b0 = rep(1, nrow(st_coordinates(area.pre[,1])))), s = indexs)
  )
  
  return(inla.stack(stk.e, stk.p))
  
}
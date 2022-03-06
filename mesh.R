
Mesh = function(bd.sf,max.edge = c(0.7, 0.7),cut.off = 0.25,offset = c(-0.05, -0.05)){
  bd.sp = as(bd.sf,Class = 'Spatial')
  inla.mesh.2d(boundary = bd.sp,
               max.edge = max.edge, cutoff = cutoff,
               offset = offset)
}

#**************Generate Prediction DataFrame************
pred = function(res,stack.full,area.pre){
  index <- inla.stack.index(stack.full, tag = "pred")$data
  
  # Prediction
  
  pred_mean <- res$summary.fitted.values[index, "mean"]
  pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  
  coop = as.matrix (st_coordinates(area.pre[,1]))
  
  pre_r = data.frame(x = coop[, 1], y = coop[, 2], pred_mean, pred_ll,pred_ul)
  pre_c <- rbind(
    data.frame(
      x = coop[, 1], y = coop[, 2],
      value = pred_mean, variable = "pred_mean"
    ),
    data.frame(
      x = coop[, 1], y = coop[, 2],
      value = pred_ll, variable = "pred_ll"
    ),
    data.frame(
      x = coop[, 1], y = coop[, 2],
      value = pred_ul, variable = "pred_ul"
    )
  )
  pre_c$variable <- as.factor(pre_c$variable)
  list(pre_r,pre_c)
}



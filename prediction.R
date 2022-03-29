pred = function(res,stack.full,target,areaPre.sf = NULL,proN){
  index <- inla.stack.index(stack.full, tag = "pred")$data
  
  # Prediction
  
  pred_mean <- res$summary.fitted.values[index, "mean"]
  pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  
  coop = as.matrix (st_coordinates(target[,1]))
  
  pre_r = data.frame(px = coop[, 1], py = coop[, 2], pred_mean, pred_ll,pred_ul)
  pre.sf = pre_r %>%
    st_as_sf(coords = c("px", "py"), dim = "XY") %>% 
    st_set_crs(proN) %>%
    st_cast("MULTIPOINT")
  
  if(is.null(areaPre.sf) == F) {
    group = areaPre.sf$geometry
    pre.sf = aggregate(pre.sf, group, mean)
  }
  return(pre.sf)
}



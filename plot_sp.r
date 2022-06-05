plot_sp <- function(results){
  coop <- st_coordinates(results)
  pred_mean = results$pred_mean
  pred_ll = results$pred_ll
  pred_ul = results$pred_ul
  
  dpm <- rbind(
    data.frame(
      X = coop[, 1], Y = coop[, 2],
      PM2.5 = pred_mean, variable = "Mean"
    ),
    data.frame(
      X = coop[, 1], Y = coop[, 2],
      PM2.5 = pred_ll, variable = "2.5% "
    ),
    data.frame(
      X = coop[, 1], Y = coop[, 2],
      PM2.5 = pred_ul, variable = "97.5%"
    )
  )
  dpm$variable <- as.factor(dpm$variable)
  
  p = ggplot(dpm, aes(X,Y,color = PM2.5)) +
    geom_tile(size = 1)+
    facet_grid(~ variable)+
    coord_fixed(ratio = 1)+
    geom_sf(data = boundaryregion, inherit.aes = F, fill = NA)
  
 return(p) 
}
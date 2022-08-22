library(blockCV)
library(caret)
library(gstat)

fnCreateCV <- function(depoint, k, n, type, block_K = F) {
  fnCheckCreateCV(depoint, type) 
  
  if(type == 'standard'){flds <- createFolds(depoint$value, k = k, list = TRUE, returnTrain = FALSE)}
  
  if(type == 'block'){
    coop <- st_coordinates(depoint)[, c(1, 2)]
    data <- as.data.frame(cbind(coop, depoint$value))
    coordinates(data) <- ~ X + Y
    TheVariogram <- variogram(V3 ~ 1, data = data)
    d <- TheVariogram[which(TheVariogram$gamma == max(TheVariogram$gamma)), ]$dist # find distance of bound of spatial autocorr
  
    ############  Generate spatial clusters to be used in spatial K-fold Cross-Validation (Spatial K-fold CV)
    mdist <- dist(coop) # distance
    # summary(mdist)
    hc <- hclust(mdist, method = "complete") # Hierarchical Clustering based on the distance
    if(!block_K){ depoint$Clust <- cutree(hc, h = n * d) # cut the clusters based on our specified distance 
    }
    
    # cut tree based on the number of blocks, similar to k fold CV
    if(block_K){depoint$Clust <- cutree(hc, k = k)}
    index <- unique(depoint$Clust)
    
    flds <- list()
    for(i in 1:length(index)){
      id <- which(depoint$Clust == levels(as.factor(index))[i])
      flds[length(flds)+1] <- list(id)
    }
  } 
    
  if(type == 'buffering'){
    bf <- buffering(
      speciesData = depoint,
      theRange = n *d
    ) 
    flds <- bf$folds
    }

  
 return(flds)
}


fnCheckCreateCV(depoint) <- function(depoint, type){
  if(class(depoint)[[1]] != "sf"){
    stop("'The input depoint should be 'sf' obj")
  }
  
  if(sum(type %in% c('standard', 'block', 'buffering')) != 1){
    stop("'Please Check the type of the CV'")
  }
  
  if(is.null(depoint) == F && sum(c("value", "geometry") %in% colnames(depoint)) != 2 ){
    stop("'depoint' must have 'value','geometry' as column names")
  }
}
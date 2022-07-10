#####################################################################
# the file is only for generating training data and validation data #
#####################################################################

library(raster)
library(RandomFields)
library(spatstat)
library(stars)
library(sf)
library(ggplot2)
library(INLA)
library(blockCV)
library(caret)
library(gstat)

source("samplegenerator.r")
crsproj <- 4326

R <- ppp(x = runif(100), y = runif(100),c(0,1),c(0,1))
Win <- Window(R)
boundaryregion = sf::st_as_sf(Win)%>%sf::st_set_crs(crsproj)
plot(boundaryregion)

#----
# Generate surface
# Region study unit square
r <- raster(nrows = 100, ncols = 100, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
x <- coordinates(r)[, 1]
y <- coordinates(r)[, 2]


model <- RMwhittle(scale = 0.05, nu = 2, var = 1)
simu <- RFsimulate(model = model, x = x, y = y)
values(r) <- simu$variable1
plot(r)

p1 <- punifsample(100, r, crsproj)

# standard cross validation
#----
k <- 5
flds <- createFolds(p1$value, k = k, list = TRUE, returnTrain = FALSE)
for (j in 1:k) {
  id <- flds[[j]]
  training_data <- p1[-id, ]
  validation_data <- p1[id, ]
  # Next step: train the model
}
#
# Example for check
id <- flds[[1]]
training_data <- p1[-id, ]
validation_data <- p1[id, ]
ggplot() +
  geom_sf(data = training_data) +
  geom_sf(data = validation_data, aes(color = "red"))
#----
# spatial block
#----
# Computes sample (empirical) variogram
coop <- st_coordinates(p1)[, c(1, 2)]
data <- as.data.frame(cbind(coop, p1$value))
coordinates(data) <- ~ X + Y
TheVariogram <- variogram(V3 ~ 1, data = data)
d <- TheVariogram[which(TheVariogram$gamma == max(TheVariogram$gamma)), ]$dist # find distance of bound of spatial autocorr

############
############  Generate spatial clusters to be used in spatial K-fold Cross-Validation (Spatial K-fold CV)
############
mdist <- dist(coop) # distance
# summary(mdist)
hc <- hclust(mdist, method = "complete") # Hierarchical Clustering based on the distance

# for spatial block CV, h is larger than d
p1$Clust <- cutree(hc, h = 1.2 * d) # cut the clusters based on our specified distance

flds <- unique(p1$Clust)

# perform CV

for (j in 1:length(flds)) {
  id <- which(p1$Clust == levels(as.factor(flds))[j])
  training_data <- p1[-id, ]
  validation_data <- p1[id, ]
}

# Example for check
id <- which(p1$Clust == levels(as.factor(flds))[1])
training_data <- p1[-id, ]
validation_data <- p1[id, ]
ggplot() +
  geom_sf(data = training_data) +
  geom_sf(data = validation_data, aes(color = "red"))
# Remark: We don't specify the how many blocks that we want, but focus more on
# really define a spatial independent cluster, if we want specify K, we can use p1$Clust = cutree(hc, k = K). Then it will give K fold, but the
# block size can differ. The blockCV's function is as follows, but when I check the results, some of them are wired. Maybe it is because it is for bio-data


sb1 <- spatialBlock(
  speciesData = p1,
  theRange = d * 111325,
  k = 5
)

training_data <- p1[sb1$folds[1][[1]][[1]], ]
validation_data <- p1[sb1$folds[1][[1]][[2]], ]

ggplot() +
  geom_sf(data = training_data) +
  geom_sf(data = validation_data, aes(color = "red"))
#----
# buffering
#----
# Here I prefer to select a buffering points because if we have many point data, e.g. 100/1000, the computational cost can be very high,
# for buffering, the range is just d
N <- 10
for (i in 1:N) {
  id_test <- sample(seq(1, nrow(p1), 1), 1)

  bf <- buffering(
    speciesData = p1,
    theRange = 0.8*d * 111325
  ) # change WSG84 to meters

  training_data <- p1[bf$folds[id_test][[1]][[1]], ]
  validation_data <- p1[bf$folds[id_test][[1]][[2]], ]
}

# check
ggplot() +
  geom_sf(data = training_data) +
  geom_sf(data = validation_data, aes(color = "red"))

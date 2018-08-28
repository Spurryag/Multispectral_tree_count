#Import libraries---------
library(raster)
library(rgdal)
library(rgeos)
library(landsat)
library(RStoolbox)
library(cluster)
library(NbClust)
library(factoextra)
library(tmap)
library(doParallel)
library(randomForest)
library(caret)
library(ROCR)
library(stargazer)
library(simpleboot)

#Set seed and multicore-------
set.seed(12345)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

#Plot the different light bands from the Orthonomic photos -------------
#below code sourced from :
#https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/landsat-data-in-r-geotiff/

#Color band combinations for landsat 5:
#https://www.esri.com/arcgis-blog/products/product/imagery/band-combinations-for-landsat-8/?rmedium=redirect&rsource=/esri/arcgis/2013/07/24/band-combinations-for-landsat-8

spec_rgb <-
  stack(
    "C:/Users/Sanil Purryag/Google Drive/St Andrews docs/Modules/Semester 2/Thesis data/PortMoak_RTK_Ortho.tif"
  )

spec_br <- brick(spec_rgb)

#rename the different bands
names(spec_rgb) <- c('Blue', 'Green', 'Red', 'NIR', 'RE')

names(spec_br) <- c('Blue', 'Green', 'Red', 'NIR', 'RE')

#check if files have been properly renamed

names(spec_br)

#Import the DSM data
spec_lidar <-
  raster(
    "C:/Users/Sanil Purryag/Google Drive/St Andrews docs/Modules/Semester 2/Thesis data/PortMoak_RTK_DSM.tif"
  )

#Plot each individual band

plot(
  spec_br,
  stretch = "hist",
  maxv = 65535,
  #col = gray(0:100 / 100),
  box = FALSE,
  axes = FALSE
)

#True color composite plot

par(col.axis = "white",
    col.lab = "white",
    tck = 0)
plotRGB(
  spec_br,
  r = 3,
  g = 2,
  #should be 2.5 times more
  b = 1,
  stretch = "hist",
  axes = FALSE,
  scale = 65535,
  main = "RGB composite image/n Landsat Bands 3,2,1"
)


#False color composite plot
par(col.axis = "white",
    col.lab = "white",
    tck = 0)
plotRGB(
  spec_br,
  r = 5,
  g = 4,
  b = 3,
  stretch = "hist",
  axes = FALSE,
  main = "RGB False composite image\n Landsat Bands 5,4,3"
)
box(col = "white")


#Exploratory analysis------------

#Plot histogram of values from DSM

#assign and replace NAs with base 160
histlidar <- spec_lidar
NAvalue(histlidar) <- 160
#Plot distribution of values in hist
hist(
  histlidar,
  main = "Distribution of DSM Values(100,000 pixels sampled)",
  xlab = "DSM Elevation Value (Inches)",
  ylab = "Frequency",
  col = "cornflowerblue",
  ylim = c(0, 20000),
  xlim = c(150, 200)
)

hist(spec_br,
     xlab = "Pixel Intensity",
     ylab = "Frequency",
     col = "cornflowerblue")
# ylim = c(0, 20000),
# xlim = c(150, 200))

#Pairs plot of different bands
pairsplot <- pairs(spec_rgb)

#Hillshading (DSM Workstream)-----------
#Below code adapted from:
#https://github.com/mgimond/Spatial/blob/master/A4-Raster-Operations.Rmd
#https://www.rdocumentation.org/packages/raster/versions/2.6-7/topics/focal

#Estimate the slope and aspect of the terrain
slope <- terrain(spec_lidar, opt = 'slope')
aspect <- terrain(spec_lidar, opt = 'aspect')

#Create hillshade (DTM) using the queen case (neighbours = 8)
#Current specification used for rough terrain estimation
hill <- hillShade(slope, aspect, 40, 270)

#plot the tmap package to see greyscale raster
hill <-
  aggregate(hill, 25) # aggregate the raster by a factor of 25
tm_shape(hill) + tm_raster(palette = "Greys") +
  tm_legend(outside = TRUE, text.size = .6) + tm_layout(title = "Hillshaded DSM")

#Histogram of hillside values - 100,000 sampled values
hist(
  hill,
  main = "Distribution of Aggregated hill values",
  xlab = "Hillshaded values",
  ylab = "Frequency",
  col = "cornflowerblue"
)
abline(v = 0.5, col = "Red", lty = 4)

#Zonal statistics - Sobel----

#===============================================================================Sobel Kernel Application

#Assign hill to new object
hillside <- hill
#Threshold hillside
hillside[hillside > 0.5] <- NA
na.omit(hillside)
#plot code adapted from:
#https://stackoverflow.com/questions/28247812/r-raster-package-plot-producing-monochrome-images

plot(hill,
     main = "Hillshaded DSM raster",
     col = gray.colors(
       10,
       start = 0.3,
       end = 0.9,
       gamma = 2.2,
       alpha = NULL
     ))
plot(
  hillside,
  main = "Thresholded DSM hillshaded raster at 0.5 level",
  col = gray.colors(
    10,
    start = 0.3,
    end = 0.9,
    gamma = 2.2,
    alpha = NULL
  )
)

#Sobel weight matrix definition
Sobel <- matrix(c(-1, 0, 1,-2, 0, 2,-1, 0, 1) / 4, nrow = 3)
f5    <- focal(hillside, w = Sobel, fun = sum)
# tm_shape(f5) + tm_raster(palette = "black") +
#   tm_legend(legend.show = TRUE, "Distribution of values")
freq(f5)
plot(f5,
     col = gray.colors(
       10,
       start = 0.3,
       end = 0.9,
       gamma = 2.2,
       alpha = NULL
     ),
     main = "Sobel Edge Detection Application")
#Flatten the map by reprojecting the raster to a Km by km dimension and to the British grid
Projected.sobel <-
  projectRaster(f5, crs = "+init=epsg:27700 +units=km")
#project the raster to a resolution of 10m per cell
Reprojection.sobel <-
  projectRaster(Projected.sobel,
                crs = crs(Projected.sobel),
                res = 0.010)
plot(Reprojection.sobel, main = "Reprojected Sobel Raster")
count.sobel <- freq(Reprojection.sobel)
print(count.sobel)

#===============================================================================Sobel over different levels
# #Different levels of thresholding based on 0.01 increase between 0.1 and 1
# #100 numbers between 0 and 1
# seq.length <- seq(0.1, 1, 0.01)
# sobel.count <- list()
#
# for (i in seq.length) {
#   #Reset before thresholding
#   hillside <- hill
#   #threshold hillside by i
#   hillside[hillside > i] <- NA
#   na.omit(hillside)
#   #Compute sobel kernel and associated focal stats
#   Sobel <- matrix(c(-1, 0, 1,-2, 0, 2,-1, 0, 1) / 4, nrow = 3)
#   f6    <- focal(hillside, w = Sobel, fun = sum)
#   #Reproject the computed focal raster
#   Projected.sobel <-
#     projectRaster(f6, crs = "+init=epsg:27700 +units=km")
#   #project the raster to a resolution of 10m per cell
#   Reprojection.sobel <-
#     projectRaster(Projected.sobel,
#                   crs = crs(Projected.sobel),
#                   res = 0.010)
#   #Obtain the frequency of observed values
#   count.sobel <- as.matrix(freq(Reprojection.sobel))
#   sobel.count[[paste0("Threshold level ", i)]] <- count.sobel
#   #print(sobel.count)
# }
#
# tree.sobel <- c()
# for (m in 1:length(sobel.count)) {
#   tree.sobel <-
#     c(tree.sobel, (sobel.count[[m]][(as.numeric(which(sobel.count[[m]] == 0))), 2]))
# }
#
# plot(
#   y = tree.sobel,
#   x = seq.length,
#   type = "o",
#   ylab = "Count of Trees",
#   xlab = "Thresholding Values",
#   main = "Sobel Tree count evolution based on Thresholding ",
#   col = "cornflowerblue"
# )
#===============================================================================Sobel Adjusted thresholds
#Different levels of thresholding based on 0.01 increase between 0 and 0.6
#100 numbers between 0 and 1
seq.length1 <- seq(0.1, 0.6, 0.001)
sobel.count5 <- list()

for (h in seq.length1) {
  #Reset before thresholding
  hillside <- hill
  #threshold hillside by i
  hillside[hillside > h] <- NA
  na.omit(hillside)
  #Compute sobel kernel and associated focal stats
  Sobel <- matrix(c(-1, 0, 1,-2, 0, 2,-1, 0, 1) / 4, nrow = 3)
  f7    <- focal(hillside, w = Sobel, fun = sum)
  #Reproject the computed focal raster
  Projected.sobel5 <-
    projectRaster(f7, crs = "+init=epsg:27700 +units=km")
  #project the raster to a resolution of 10m per cell
  Reprojection.sobel5 <-
    projectRaster(Projected.sobel5,
                  crs = crs(Projected.sobel5),
                  res = 0.010)
  #Obtain the frequency of observed values
  count.sobel5 <- as.matrix(freq(Reprojection.sobel5))
  sobel.count5[[paste0("Threshold level ", h)]] <- count.sobel5
  #print(sobel.count5)
}

tree.sobel1 <- c()
for (o in 1:length(sobel.count5)) {
  tree.sobel1 <-
    c(tree.sobel1, (sobel.count5[[o]][(as.numeric(which(sobel.count5[[o]] == 0))), 2]))
}

plot(
  y = tree.sobel1,
  x = seq.length1,
  type = "o",
  ylab = "Count of Trees",
  xlab = "Thresholding Values",
  main = "Sobel Tree count evolution based on Thresholding ",
  col = "cornflowerblue"
)



#===============================================================================Sobel Non parametric CI

#compute the SE of the simulated distribution
standard.error.sobel <- sd(tree.sobel1) / sqrt(length(tree.sobel1))
mean.error.sobel  <- mean(tree.sobel1)
#upper and lower normal limits
ul.sobel <- mean.error.sobel + 1.96 * (standard.error.sobel)
ll.sobel <- mean.error.sobel -  1.96 * (standard.error.sobel)
#add lines for bands to the plot y axis
abline(h = mean.error.sobel, col = "black", lty = 2)
abline(h = ul.sobel, col = "red", lty = 2)
abline(h = ll.sobel, col = "red", lty = 2)

#compute a nonparametric bootstrap from the different thresholded counts
#code adapted from:
#https://stats.stackexchange.com/questions/16516/how-can-i-calculate-the-confidence-interval-of-a-mean-in-a-non-normally-distribu

#Take tree counts from different thresholds and laplace application
#2000 resamples, 5% trim

sobel.boot <- one.boot(tree.sobel1, mean, R = 2000, tr = 0.05)
hist(sobel.boot, main = "Sobel count bootstrap distribution", col = "cornflowerblue")
sobel.bootstrap <- boot.ci(sobel.boot, type = c("perc", "bca"))
print(sobel.bootstrap)


#Zonal statistics - Laplace----

#===============================================================================Laplacian filter Application

#Assign hill to new object
hillside <- hill
#Threshold hillside
hillside[hillside > 0.5] <- NA
na.omit(hillside)

laplace <- matrix(c(0, 1, 0, 1,-4, 1, 0, 1, 0), nrow = 3)
f10    <- focal(hillside, w = laplace, fun = sum)
#tm_shape(f10) + tm_raster(palette = "black") +
#  tm_legend(legend.show = TRUE)
freq(f10)
plot(f10,
     col = gray.colors(
       10,
       start = 0.3,
       end = 0.9,
       gamma = 2.2,
       alpha = NULL
     ),
     main = "Laplace Edge Detection Application")

#Flatten the map by reprojecting the raster to a Km by km dimension and to the British grid
Projected.laplace <-
  projectRaster(f10, crs = "+init=epsg:27700 +units=km")
#project the raster to a resolution of 10m per cell
Reprojection.laplace <-
  projectRaster(Projected.laplace,
                crs = crs(Projected.laplace),
                res = 0.010)
count.laplace <- freq(Reprojection.laplace)
print(count.laplace)
plot(Reprojection.laplace, main = "Reprojected Laplace")
#===============================================================================Laplacian over different levels

# #Different levels of thresholding based on 0.01 increase between 0.1 and 1
# #100 numbers between 0 and 1
# seq.length <- seq(0.1, 1, 0.01)
# laplace.count <- list()
#
# for (i in seq.length) {
#   #Reset before thresholding
#   hillside <- hill
#   #threshold hillside by i
#   hillside[hillside > i] <- NA
#   na.omit(hillside)
#   #Compute sobel kernel and associated focal stats
#   laplace <- matrix(c(0, 1, 0, 1,-4, 1, 0, 1, 0), nrow = 3)
#   f11    <- focal(hillside, w = laplace, fun = sum)
#   #Reproject the computed focal raster
#   Projected.laplace <-
#     projectRaster(f11, crs = "+init=epsg:27700 +units=km")
#   #project the raster to a resolution of 10m per cell
#   Reprojection.laplace <-
#     projectRaster(Projected.laplace,
#                   crs = crs(Projected.laplace),
#                   res = 0.010)
#   #Obtain the frequency of observed values
#   count.laplace <- as.matrix(freq(Reprojection.laplace))
#   laplace.count[[paste0("Threshold level ", i)]] <- count.laplace
#   #print(laplace.count)
# }
#
# tree.laplace <- c()
# for (m in 1:length(laplace.count)) {
#   tree.laplace <-
#     c(tree.laplace, (laplace.count[[m]][(as.numeric(which(laplace.count[[m]] ==
#                                                             0))), 2]))
# }
#
# plot(
#   y = tree.laplace,
#   x = seq.length,
#   type = "o",
#   ylab = "Count of Trees",
#   xlab = "Thresholding Values",
#   main = "Laplace Tree count evolution based on Thresholding ",
#   col = "cornflowerblue"
# )
#===============================================================================Laplacian over bounded threshold
#Different levels of thresholding based on 0.01 increase between 0.1 and 0.6
#100 numbers between 0 and 1
seq.length2 <- seq(0.1, 0.6, 0.001)
laplace.count5 <- list()

for (y in seq.length2) {
  #Reset before thresholding
  hillside <- hill
  #threshold hillside by i
  hillside[hillside > y] <- NA
  na.omit(hillside)
  #Compute sobel kernel and associated focal stats
  laplace <- matrix(c(0, 1, 0, 1,-4, 1, 0, 1, 0), nrow = 3)
  f16    <- focal(hillside, w = laplace, fun = sum)
  #Reproject the computed focal raster
  Projected.laplace5 <-
    projectRaster(f16, crs = "+init=epsg:27700 +units=km")
  #project the raster to a resolution of 10m per cell
  Reprojection.laplace5 <-
    projectRaster(Projected.laplace5,
                  crs = crs(Projected.laplace5),
                  res = 0.010)
  #Obtain the frequency of observed values
  count.laplace5 <- as.matrix(freq(Reprojection.laplace5))
  laplace.count5[[paste0("Threshold level ", y)]] <- count.laplace5
  #print(laplace.count)
}

tree.laplace1 <- c()
for (d in 1:length(laplace.count5)) {
  tree.laplace1 <-
    c(tree.laplace1, (laplace.count5[[d]][(as.numeric(which(laplace.count5[[d]] ==
                                                              0))), 2]))
}

plot(
  y = tree.laplace1,
  x = seq.length2,
  type = "o",
  ylab = "Count of Trees",
  xlab = "Thresholding Values",
  main = "Laplace Tree count evolution based on Thresholding ",
  col = "cornflowerblue"
)

#compute the SE of the simulated distribution
standard.error.laplace <-
  sd(tree.laplace1) / sqrt(length(tree.laplace1))
mean.error.laplace <- mean(tree.laplace1)
#upper and lower normal limits
ul.laplace <- mean.error.laplace + 1.96 * (standard.error.laplace)
ll.laplace <- mean.error.laplace - 1.96 * (standard.error.laplace)
#add lines for bands to the plot y axis
abline(h = mean.error.laplace, col = "black", lty = 1)
abline(h = ul.laplace, col = "red", lty = 2)
abline(h = ll.laplace, col = "red", lty = 2)


#compute a nonparametric bootstrap from the different thresholded counts
#code adapted from:
#https://stats.stackexchange.com/questions/16516/how-can-i-calculate-the-confidence-interval-of-a-mean-in-a-non-normally-distribu

#Take tree counts from different thresholds and laplace application
#2000 resamples, 5% trim

laplace.boot <- one.boot(tree.laplace1, mean, R = 2000, tr = 0.05)
hist(laplace.boot, main = "laplace count bootstrap distribution", col = "cornflowerblue")
laplace.bootstrap <- boot.ci(laplace.boot, type = c("perc", "bca"))
print(laplace.bootstrap)

#Measures of accuracy for zonal stats---------
#Code adapted from:
#https://gis.stackexchange.com/questions/119305/how-to-calculate-and-compare-rmse-between-two-dems
#https://gis.stackexchange.com/questions/4802/rmse-between-two-rasters-step-by-step

#Extract the values from the thresholded rasters at 0.5
hill.error <- getValues(hillside)
f5.error <- getValues(f5)
f10.error <- getValues(f10)

#Sobel at 0.5 threshold (RMSE)
#sqrt(mean(Actual - Predicted)^2)
diff.sobel <- (hill.error - f5.error) ^ 2
avg.sobel <- mean(na.omit(diff.sobel))
#RMSE.sobel <- sqrt(avg.sobel)
print(avg.sobel)

#Sobel at 0.5 threshold (PNSR)
PNSR.sobel <- 10 * log(((0.6910431 ^ 2) / avg.sobel))
print(PNSR.sobel)

#Laplacian at 0.5 threshold (RMSE)
diff.laplace <- (hill.error - f10.error) ^ 2
avg.laplace <- mean(na.omit(diff.laplace))
#RMSE.laplace <- sqrt(avg.laplace)
print(avg.laplace)

#Laplacian at 0.5 Threshold (PNSR)
PNSR.laplace <- 10 * log(((2.675293 ^ 2) / avg.laplace))
print(PNSR.laplace)

#Compute the different index and color composites--------------

#Aggregate Orthophoto by a factor of 45
spec_br1 <- aggregate(spec_rgb, 45, progress = 'text')

names(spec_br1) <- c('Blue', 'Green', 'Red', 'NIR', 'RE')

#Normalized Differentiation Vegetation Index (NDVI)
#Formula sourced from:
##https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/vegetation-indices-NDVI-in-R/

NIR <- spec_br1[[4]]
Red <- spec_br1[[3]]
Blue <- spec_br1[[1]]

#Compute the NDVI raster
Spec_ndvi <-
  overlay(
    NIR,
    Red,
    fun = function(x, y) {
      (x - y) / (x + y)
    }
  )

#Plot the NDVI
par(mar = c(1, 1, 1, 0.5))
plot(Spec_ndvi, main = "Normalized differentiation vegetation index", box = FALSE)

#Enhanced Vegetation index (EVI)
#Formula sourced from:
#https://cran.r-project.org/web/packages/LSRS/LSRS.pdf

Spec_evi <- overlay(
  NIR,
  Red,
  Blue,
  fun = function(x, y, z) {
    2.5 * ((x - y) / (x + (6 * y) + (7.5 * z) + 1))
  }
)

#Unsupervised machine learning on stacked orthophoto (Kmeans)---------------
#Formula:
#https://medium.com/regen-network/remote-sensing-indices-389153e3d947

#Below code adapted from:
#http://rspatial.org/analysis/rst/9-remotesensing.html
#https://geoscripting-wur.github.io/AdvancedRasterAnalysis/


#Create determinant rasters of the cluster
#benchmark <- stack(Spec_evi, Spec_ndvi, spec_br1[[5]])
benchmark <- stack(Spec_evi, Spec_ndvi)

#==============================================================================Plot all clustering possibilities
kposs <- c()
kpossibilities <- function(num = 10) {
  for (k in 1:num) {
    set.seed(12345)
    kmncluster <-
      kmeans(
        na.omit(benchmark[]),
        centers = k,
        iter.max = 500,
        nstart = 5,
        algorithm = "Lloyd"
      )
    kposs <- c(kposs, kmncluster$tot.withinss)
    kclass <- spec_br1[[1]]
    kclass[] <- kmncluster$cluster
    plot(kclass, main = paste(k, "cluster(s)"))
  }
  return(kposs)
}

totalwss <- kpossibilities()
plot(
  y = totalwss,
  x = seq(1, 10, 1),
  main = "Elbow Method Results",
  type = "o",
  ylab = "Total Within Sum of Squares",
  xlab  = "Number of clusters",
  col = "cornflowerblue"
)
abline(v = 3,
       lty = 2,
       col = "red",
       lwd = 3)
abline(h = 1840.6808,
       lty = 2,
       col = "red",
       lwd = 3)


# nbclustering.avg <-
#   NbClust(
#     data = benchmark[],
#     diss = NULL,
#     distance = "euclidean",
#     min.nc = 2,
#     max.nc = 15,
#     method = "average",
#     index = "cindex",
#     alphaBeale = 0.1
#   )

#============================================================================== Plot desired and create the training data

#Attempt to perform 3 cluster classification based on reduced NDVI image

#Kmeans clustering result
set.seed(12345)
kmncluster <-
  kmeans(
    na.omit(benchmark[]),
    centers = 3,
    iter.max = 500,
    nstart = 5,
    algorithm = "Lloyd"
  )

#Create training polygons and assign to arbitrary raster for visualisation
kclass <- spec_br1[[5]]
kclass[] <- kmncluster$cluster
#plot clustering result
plot(kclass)
freq(kclass)

#Edit the clustered raster details
classes <- kclass
names(classes) <- "Classes"
#plot Thematic raster
par(mar = c(1, 1, 1, 0.5))
plot(classes, main = "Thematic Raster")

#Create training data by isolating unwanted regions/classes through thresholding
#1 and 2 appears to be redundant here
classification <- classes
#plot(classification)
#Clusters less than 6 are omitted
#classification[classification < 6] <- NA
classification[classification > 2] <- NA

plot(classification)

#isolate the predictive bands; bands 1-4
#preds <- stack(spec_br1[[1]],spec_br1[[2]],spec_br1[[3]],spec_br1[[4]])
preds <- spec_br1
#mask predictors with defined classification to obtain training data
covmasked <- mask(preds, classification)
plot(covmasked)
trainlayer <- addLayer(covmasked, classification)
plot(trainlayer)

#Extract the values from training stack
values <- getValues(trainlayer)
valuetable <- na.omit(values)
valuetable <- as.data.frame(valuetable)
#colnames(valuetable)[5] <- "Classes"
valuetable$Classes <- as.factor(valuetable$Classes)
print(unique(valuetable$Classes))

#All band possibilities#####

allband <- stack(spec_br1, Spec_ndvi, Spec_evi)
names(allband) <-
  c('Blue', 'Green', 'Red', 'NIR', 'RE', 'NDVI', 'EVI')


kpossibilities1 <- function(num = 10) {
  for (k in 1:num) {
    set.seed(12345)
    kmncluster <-
      kmeans(
        na.omit(allband[]),
        centers = k,
        iter.max = 500,
        nstart = 5,
        algorithm = "Lloyd"
      )
    
    kclass <- spec_br1[[1]]
    kclass[] <- kmncluster$cluster
    plot(kclass, main = paste(k, "cluster(s)"))
  }
}

kpossibilities1()


#Graph the valuetable based on the associated classes####
#Plot the different bands across the classes####
band.clust <- function(cluster) {
  par(mfrow = c(2, 3))
  val_class <- subset(valuetable, Classes == cluster)
  for (p in 1:5) {
    hist(
      val_class[[p]],
      main = c(paste(
        "Cluster",
        paste(cluster),
        colnames(val_class)[p],
        " Band distrbution"
      )),
      col = "cornflowerblue",
      xlab = c(paste(colnames(val_class)[p]), " values")
    )
  }
}


for (v in 6:7) {
  band.clust(v)
}

dev.off()

#Compute a tree count using kmeans-------------

classification <- classes

Projected.kmeans <-
  projectRaster(classification, crs = "+init=epsg:27700 +units=km", res = 0.010)

plot(Projected.kmeans, main = "kmeans projected raster")

freq(Projected.kmeans)

#Supervised machine learning on stacked orthophoto (RandomForest)---------------
#below code adapted from:
##https://geoscripting-wur.github.io/AdvancedRasterAnalysis/

tune <- tuneRF(x = valuetable[, c(1:5)],
               y = valuetable$Classes)
#plot OOB error
plot(
  tune,
  type = "o",
  main = "RF Tuning result",
  ylab = "OOB Error percentage",
  xlab = "Number of variables randomly sampled at each split",
  col = "cornflowerblue"
)

set.seed(12345)
#Fit the Random forest model
modelRF.final <-
  randomForest(
    x = valuetable[, c(1:5)],
    y = valuetable$Classes,
    mtry = 4,
    importance = TRUE
  )

#variable importance plot
varimp <- varImpPlot(modelRF.final)

#Predict tree cover on predictors (insample prediction)
predTree <- predict(preds, model = modelRF.final)

#plot predicted raster
par(mar = c(1, 1, 1, 0.5))
plot(predTree, main = "Random Forest Predicted Raster")

#Plot classes is one for tree and zero for notree
#Reproject ther raster
Projected.pred <-
  projectRaster(predTree, crs = "+init=epsg:27700 +units=km")

#project the raster to a resolution of 10m per cell
Reprojection.pred <-
  projectRaster(Projected.pred,
                crs = crs(Projected.pred),
                res = 0.010)
plot(Reprojection.pred)
#Threshold the reprojection for the count of trees (Based on color from legend)
Thresholded <- Reprojection.pred
RFcount <- freq(Thresholded)
print(RFcount)

#Edit confusion matrix for better interpretation - from computed RF model
colnames(modelRF.final$confusion) <-
  c("Trees", "Non-Trees", "class Error")
rownames(modelRF.final$confusion) <- c("Trees", "Non-Trees")
modelRF.final$confusion

#Plot sample tree from random forest
#(getTree(modelRF.final,4, labelVar = TRUE))


#below code adapted from:
#https://stackoverflow.com/questions/30366143/how-to-compute-roc-and-auc-under-roc-after-training-using-caret-in-r

#Compute the accuracy metrics
predictions <-
  as.vector(modelRF.final$votes[, 2]) # Random forest output
predictions.auc <- prediction(predictions, valuetable$Classes)
#AUC
auc.performance <- performance(predictions.auc, "auc")
#ROC
roc.performance <- performance(predictions.auc, "tpr", "fpr")
#plot the ROC
#below code adapted from:
#https://chandramanitiwary.wordpress.com/2014/03/17/r-tips-part2-rocr-example-with-randomforest/

plot(roc.performance,
     main = "ROC curve for Fitted Random Forest",
     col = "cornflowerblue",
     lwd = 2,
     xlab = "Specificity",
     ylab = "Sensitivity")
abline(
  a = 0,
  b = 1,
  lwd = 2,
  lty = 2,
  col = "gray"
)

#AUC detailed stats
auc <- unlist(slot(auc.performance, "y.values"))
#Compute min and max values
minauc <- min(round(auc, digits = 2))
maxauc <- max(round(auc, digits = 2))
minauct <- paste(c("min(AUC) = "), minauc, sep = "")
maxauct <- paste(c("max(AUC) = "), maxauc, sep = "")


#Threshold definition adapted from:
#https://stackoverflow.com/questions/16347507/obtaining-threshold-values-from-a-roc-curve
cutoffs <- data.frame(cut=roc.performance @alpha.values[[1]], fpr=roc.performance @x.values[[1]], 
                      tpr=roc.performance @y.values[[1]])
head(cutoffs)
#order cutoffs
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]

#####################################################################10% threshold

set.seed(12345)
#Fit the Random forest model
modelRF.final <-
  randomForest(
    x = valuetable[, c(1:5)],
    y = valuetable$Classes,
    mtry = 4,
    importance = TRUE,
    cutoff = c(0.9,0.1)
  )

#variable importance plot
varimp <- varImpPlot(modelRF.final)

#Predict tree cover on predictors (insample prediction)
predTree <- predict(preds, model = modelRF.final)

#plot predicted raster
par(mar = c(1, 1, 1, 0.5))
plot(predTree, main = "Random Forest Predicted Raster")

#Plot classes is one for tree and zero for notree
#Reproject ther raster
Projected.pred <-
  projectRaster(predTree, crs = "+init=epsg:27700 +units=km")

#project the raster to a resolution of 10m per cell
Reprojection.pred <-
  projectRaster(Projected.pred,
                crs = crs(Projected.pred),
                res = 0.010)
plot(Reprojection.pred)
#Threshold the reprojection for the count of trees (Based on color from legend)
Thresholded <- Reprojection.pred
RFcount <- freq(Thresholded)
print(RFcount)

#Edit confusion matrix for better interpretation - from computed RF model
colnames(modelRF.final$confusion) <-
  c("Trees", "Non-Trees", "class Error")
rownames(modelRF.final$confusion) <- c("Trees", "Non-Trees")
modelRF.final$confusion

#Plot sample tree from random forest
#(getTree(modelRF.final,4, labelVar = TRUE))


#below code adapted from:
#https://stackoverflow.com/questions/30366143/how-to-compute-roc-and-auc-under-roc-after-training-using-caret-in-r

#Compute the accuracy metrics
predictions <-
  as.vector(modelRF.final$votes[, 2]) # Random forest output
predictions.auc <- prediction(predictions, valuetable$Classes)
#AUC
auc.performance <- performance(predictions.auc, "auc")
#ROC
roc.performance <- performance(predictions.auc, "tpr", "fpr")
#plot the ROC
#below code adapted from:
#https://chandramanitiwary.wordpress.com/2014/03/17/r-tips-part2-rocr-example-with-randomforest/

plot(roc.performance,
     main = "ROC curve for Fitted Random Forest",
     col = "cornflowerblue",
     lwd = 2,
     xlab = "Specificity",
     ylab = "Sensitivity")
abline(
  a = 0,
  b = 1,
  lwd = 2,
  lty = 2,
  col = "gray"
)

#AUC detailed stats
auc <- unlist(slot(auc.performance, "y.values"))
#Compute min and max values
minauc <- min(round(auc, digits = 2))
maxauc <- max(round(auc, digits = 2))
minauct <- paste(c("min(AUC) = "), minauc, sep = "")
maxauct <- paste(c("max(AUC) = "), maxauc, sep = "")


#Threshold definition adapted from:
#https://stackoverflow.com/questions/16347507/obtaining-threshold-values-from-a-roc-curve
cutoffs <- data.frame(cut=roc.performance @alpha.values[[1]], fpr=roc.performance @x.values[[1]], 
                      tpr=roc.performance @y.values[[1]])
head(cutoffs)
#order cutoffs
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]

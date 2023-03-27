#######################################################################################
#########################  DATA PREPARATION FUNCTION  #################################
#######################################################################################

#######################################################################################
#######################################################################################
###  Written by André S.J. van Proosdij & Niels Raes, 2016
###  Adapted by Sampaio, A. C. P. (1) & Cavalcante, A. de M. B. (2), 2021
###  1 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
###  2 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
#######################################################################################
#######################################################################################

#######################################################################################
# This script is used to prepare data files with spatial environmental data. Separate
# sections are written for climatic, and altitudinal data sets. A final section
# deals with analysis of multicollinearity, selection of variables and the preparation
# of PCA axes based on selected variables.
#######################################################################################

#######################################################################################
#########################  Index  #####################################################
#######################################################################################
# 1. Load packages.
# 2. Climate data from WORLDCLIM.
# 3. Altitude data at 90 m spatial resolution.
# 4. Analysis of multicollinearity & selection of variables.
# 5. Preparing PCA axes as input variables for analysis. 
#######################################################################################

#######################################################################################
#########################  1. Load packages  ##########################################
#######################################################################################

rm(list = ls(all = TRUE))
setwd("C:/Program Files/R/R-3.5.3")
getwd()
library(raster) # stack(), scale(), crop(), writeRaster(), raster(), aggregate(), raster(), mask(), and mosaic() functions
library(rgdal)
library(dismo)
library(SDMTools) # asc2dataframe() function
library(sp)
library(adehabitatHS)
library(plyr) # joining df
library(scales) # rescale() function
require(maptools) # readShapeSpatial() function
require(rgeos) # gBuffer() function
library(Hmisc) # rcorr() function
library(ade4) # dudi.pca() function
data(wrld_simpl)

# Define the extent of the study area: lon 45W 34W; lat 2S 17S
ext.CAfr <- extent(-45, -34, -17, -2)

#######################################################################################
#########################  2. Climate data from WORLDCLIM  ############################
#######################################################################################

# Download climate data from the WORLDCLIM site: www.worldclim.org (Hijmans & al., 2005).

setwd ("C:/Mnspecies/Worldclim/30sec")
getwd()

# Read each .bil file and crop it to the spatial extent of the study area.
Worldclim.files <- list.files("C:/Mnspecies/Worldclim/30sec", pattern = "[.]bil", full.names = TRUE)
for(i in Worldclim.files) {
  Worldclim.CAfr <- crop((raster(i)), ext.CAfr)
  exportname <- strsplit(basename(i), ".bil")
  # Export as .asc file, first set the directory, otherwise writeRaster will not work.
  setwd ("C:/Mnspecies/Worldclim/30sec/cutlayers")
  writeRaster(Worldclim.CAfr,
              filename  = paste(c("C:/Mnspecies/Worldclim/30sec/cutlayers/", exportname, "CAfr"), collapse = ""),
              format    = "ascii",
              NAflag    = -9999,
              overwrite = TRUE)
}

# Aggregate each .asc file to the required spatial resolution:
# 30 arcsec to 5 arcmin -> factor 10, 30 arcsec to 2.5 arcmin -> factor 5, 30 arcsec to 1 arcmin -> factor 2.
Wclim.files <- list.files("C:/Mnspecies/Worldclim/30sec/cutlayers", pattern = "[.]asc", full.names = FALSE)
# Aggregate each file and export it as .asc file.
for(i in Wclim.files) {
  setwd("C:/Mnspecies/Worldclim/30sec/cutlayers")
  raster <- raster(i)
  Wclim.1arcmin <- aggregate(raster, fact = 2, fun = mean, expand = FALSE, na.rm = TRUE, overwrite = TRUE)
  exportname <- strsplit(basename(i), ".asc")
  # First set the directory, otherwise writeRaster will not work.
  setwd("C:/Mnspecies/Worldclim/1min/aggregatelayers")
  writeRaster(Wclim.1arcmin,
              filename  = paste("C:/Mnspecies/Worldclim/1min/aggregatelayers/", exportname, "_1min", sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = TRUE)
}

#######################################################################################
#########################  3. Altitude data at 3 arcsec spatial resolution  ###########
#######################################################################################

# Download altitude data at 3 arcsec spatial resolution (approx. 90 m at the equator)
# from the SRTM site: http://srtm.csi.cgiar.org/srtmdata/

setwd ("C:/Mnspecies/DEM")
getwd()

DEM.files <- list.files("C:/Mnspecies/DEM", pattern = "[.]asc", full.names = FALSE)

#######################################################################################
# Defined altitude data.
#######################################################################################
# Based on the data at 3 arcsec spatial resolution, defined altitude data at a coarser spatial 
# resolution. Aggregate each .asc file to the required spatial resolution: 3 arcsec to 
# 5 arcmin -> factor 100, 3 arcsec to 2.5 arcmin -> factor 50, 3 arcsec to 1 arcmin -> 
# factor 20, 3 arcsec to 30 arcsec -> factor 10.

for(i in DEM.files)  {
  setwd("C:/Mnspecies/DEM")
  raster <- raster(i)
  Alt.1arcmin <- aggregate(raster, fact = 20, fun = mean, expand = FALSE, na.rm = TRUE, overwrite = TRUE)
  exportname <- strsplit(basename(i), ".asc")
  setwd("C:/Mnspecies/DEM/1min/alt") # Set the directory, otherwise writeRaster will not work
  writeRaster(Alt.1arcmin,
              filename  = paste("C:/Mnspecies/DEM/1min/alt/", exportname, "_1min", sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = TRUE)
}

# Merge / mosaic the files. Take 1 file and mosaic each following file to it.
setwd("C:/Mnspecies/DEM/1min/alt")
Alt.1min.files <- list.files("C:/Mnspecies/DEM/1min/alt", pattern = "[.]asc", full.names = FALSE)
Alt.1min.files2 <- Alt.1min.files[-1]
x <- raster(Alt.1min.files[1])
for(i in Alt.1min.files2) {
  y <- raster(i)
  z <- mosaic(x, y, fun = mean)
  x <- z
}
setwd("C:/Mnspecies/DEM/1min/alt")
writeRaster(x,
            filename  = paste("C:/Mnspecies/DEM/1min/Alt_1min_mosaic"),
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)
# Crop the mosaiced file to the extent of the study area and export as .asc file.
Alt.CAfr <- crop(x, ext.CAfr)
writeRaster(Alt.CAfr,
            filename  = "C:/Mnspecies/DEM/1min/altCAfr_1min",
            format    = "ascii",
            NAflag    = -9999,
            overwrite = TRUE)


#######################################################################################
# Calculate the standard deviation of altitude.
#######################################################################################
# Based on the data at 3 arcsec spatial resolution, calculate the standard deviation of
# altitude at a coarser spatial resolution. # Aggregate each .asc file to the required
# spatial resolution: 3 arcsec to 5 arcmin -> factor 100, 3 arcsec to 2.5 arcmin ->
# factor 50, 3 arcsec to 1 arcmin -> factor 20, 3 arcsec to 30 arcsec -> factor 10.

for(i in DEM.files)  {
  setwd("C:/Mnspecies/DEM")
  raster <- raster(i)
  DEM.1arcmin <- aggregate(raster, fact = 20, fun = sd, expand = FALSE, na.rm = TRUE, overwrite = TRUE)
  exportname <- strsplit(basename(i), ".asc")
  setwd("C:/Mnspecies/DEM/1min/sd") # Set the directory, otherwise writeRaster will not work
  writeRaster(DEM.1arcmin,
              filename  = paste("C:/Mnspecies/DEM/1min/sd/", exportname, "_1min", sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = TRUE)
}

# Merge / mosaic the files. Take 1 file and mosaic each following file to it.
setwd("C:/Mnspecies/DEM/1min/sd")
DEM.1min.sd.files <- list.files("C:/Mnspecies/DEM/1min/sd", pattern = "[.]asc", full.names = FALSE)
DEM.1min.sd.files2 <- DEM.1min.sd.files[-1]
x <- raster(DEM.1min.sd.files[1])
for(i in DEM.1min.sd.files2) {
  y <- raster(i)
  z <- mosaic(x, y, fun = mean)
  x <- z
}
setwd("C:/Mnspecies/DEM/1min/sd")
writeRaster(x,
            filename  = paste("C:/Mnspecies/DEM/1min/DEM_sd_1min_mosaic"),
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)
# Crop the mosaiced file to the extent of the study area and export as .asc file.
DEM.sd.CAfr <- crop(x, ext.CAfr)
writeRaster(DEM.sd.CAfr,
            filename  = "C:/Mnspecies/DEM/1min/DEM_sd_1min_CAfr",
            format    = "ascii",
            NAflag    = -9999,
            overwrite = TRUE)

#######################################################################################
# Calculate the range of altitude.
#######################################################################################
# Based on the data at 3 arcsec spatial resolution, calculate the range of altitude
# (maximum-minimum) at a coarser spatial resolution. # Aggregate each .asc file to the 
# required spatial resolution: 3 arcsec to 5 arcmin -> factor 100, 3 arcsec to 2.5
# arcmin -> factor 50, 3 arcsec to 1 arcmin -> factor 20, 3 arcsec to 30 arcsec -> 
# factor 10.
for(i in DEM.files)  {
  setwd("C:/Mnspecies/DEM")
  raster <- raster(i)
  DEMmax.1arcmin <- aggregate(raster, fact = 20, fun = max, expand = FALSE, na.rm = TRUE, overwrite = TRUE)
  # Aggregate the raster
  setwd("C:/Mnspecies/DEM/1min/max") # Set the directory, otherwise writeRaster will not work
  writeRaster(DEMmax.1arcmin,
              filename  = paste("C:/Mnspecies/DEM/1min/max/", i, "_1min", sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = TRUE)
  DEMmin.1arcmin <- aggregate(raster, fact = 20, fun = min, expand = FALSE, na.rm = TRUE, overwrite = TRUE)
  # Aggregate the raster
  setwd("C:/Mnspecies/DEM/1min/min")
  writeRaster(DEMmin.1arcmin,
              filename  = paste("C:/Mnspecies/DEM/1min/min/", i, "_1min", sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = TRUE)
  # Calculate the range: max-min.
  DEMmaxmin.1arcmin <- DEMmax.1arcmin - DEMmin.1arcmin
  setwd("C:/Mnspecies/DEM/1min/maxmin") # Set the directory, otherwise writeRaster will not work
  writeRaster(DEMmaxmin.1arcmin,
              filename  = paste("C:/Mnspecies/DEM/1min/maxmin/", i, "_1min", sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = TRUE)
}

# Merge / mosaic the files. Take 1 file and mosaic each following file to it.
setwd("C:/Mnspecies/DEM/1min/maxmin")
DEM.maxmin.1min.files <- list.files("C:/Mnspecies/DEM/1min/maxmin", pattern="[.]asc", full.names = FALSE)
DEM.maxmin.1min.files2 <- DEM.maxmin.1min.files[-1]
x <- raster(DEM.maxmin.1min.files[1])
for(i in DEM.maxmin.1min.files2) {
  y <- raster(i)
  z <- mosaic(x, y, fun = mean)
  x <- z
}
setwd("C:/Mnspecies/DEM/1min") # Set the directory, otherwise writeRaster will not work
writeRaster(x,
            filename  = paste("C:/Mnspecies/DEM/1min/DEM_maxmin_1min_mosaic"),
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)
# Crop the mosaiced file to the extent of the study area and export as .asc file.
DEM.max.CAfr <- crop(x, ext.CAfr)
writeRaster(DEM.max.CAfr,
            filename  = "C:/Mnspecies/DEM/1min/DEM_maxmin_1min_CAfr",
            format    = "ascii",
            NAflag    = -9999,
            overwrite = TRUE)

#######################################################################################
#########################  4. Analysis of multicollinearity & selection of variables ##
#######################################################################################

# Environmental variables are often highly correlated. To prevent errors caused by this
# multicollinearity, all variables are tested on collinearity using Spearman rank
# correlation test. From each group of correlated predictors, 1 predictor is selected.

setwd ("C:/Mnspecies/Files")
getwd()

#######################################################################################
# Calculate Spearman rank correlation for all environmental variables.
#######################################################################################
# Load all environmental variables (predictors) as .asc files. Set the directory, 
# otherwise writeRaster will not work. Put the predictors files in the defined directory
files.CAfr <- list.files('C:/Mnspecies/Files/Cropped', pattern = '.asc', full.names = TRUE)
files.CAfr.stack <- stack(files.CAfr)
files.CAfr.names <- names(files.CAfr.stack) # Get the predictor names

# Convert to dataframe, NA's are omitted, then to matrix, required by fucntion rcorr().
files.CAfr.df <- asc2dataframe(files.CAfr, varnames = files.CAfr.names)
files.CAfr.matrix <- as.matrix(files.CAfr.df)
files.CAfr.matrix <- files.CAfr.matrix[,-(1:2)] # Remove x and y columns
colnames(files.CAfr.matrix)

# Calculate the Spearman rank correlation between all variables.Set the directory, 
# otherwise writeRaster will not work.
setwd ("C:/Mnspecies/Files/Multicollinearity")
all.rcorr <- rcorr(files.CAfr.matrix, type = "spearman")
# Export the Spearman r values and significance levels.
write.table(all.rcorr$r,
            file = "all_Spearman_r.txt",
            sep = ",",
            quote = FALSE,
            append = FALSE,
            na = "NA",
            qmethod = "escape")
write.table(all.rcorr$P,
            file = "all_Spearman_significance.txt",
            sep = ",",
            quote = FALSE,
            append = FALSE,
            na = "NA",
            qmethod = "escape")

# Outside R, identify groups of correlated predictors based on Spearman rank > 0.7.
# Groups of correlated environmental variables are assessed usign the code below.

#######################################################################################
# Calculate PCA for groups of correlated environmental variables.
#######################################################################################
varnames <- colnames(files.CAfr.matrix)
varnames

# Define group 1 of collinear variables.
keep.preds1 <- c("altCAfr_1min", "bio10CAfr_1min", "bio11CAfr_1min", "bio1CAfr_1min", "bio5CAfr_1min", "bio6CAfr_1min", "bio8CAfr_1min", "bio9CAfr_1min", "bio18CAfr_1min")
files.CAfr.preds1 <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.preds1)]
# Standardize data, required for PCA analysis
files.CAfr.preds1.scale <- scale(files.CAfr.preds1, center = TRUE, scale = TRUE)
# Calculate PCA for groups of correlated predictors.
pc.preds1 <- dudi.pca(files.CAfr.preds1.scale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pc.preds1$co # the column coordinates, which is the same as vector load for each PC
write.csv(pc.preds1$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_preds1_co_CAfr.csv")
pc.preds1$eig # the eigenvalues of the PC's
barplot(pc.preds1$eig)
var1 <- pc.preds1$eig[1]/sum(pc.preds1$eig)
var1

# Define group 2 of collinear variables.
keep.preds2 <- c("bio4CAfr_1min", "bio3CAfr_1min", "bio1CAfr_1min", "bio11CAfr_1min", "bio9CAfr_1min")
files.CAfr.preds2 <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.preds2)]
# Standardize data, required for PCA analysis
files.CAfr.preds2.scale <- scale(files.CAfr.preds2, center = TRUE, scale = TRUE)
# Calculate PCA for groups of correlated predictors.
pc.preds2 <- dudi.pca(files.CAfr.preds2.scale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pc.preds2$co # the column coordinates, which is the same as vector load for each PC
write.csv(pc.preds2$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_preds2_co_CAfr.csv")

# Define group 3 of collinear variables.
keep.preds3 <- c("bio12CAfr_1min", "bio13CAfr_1min", "bio16CAfr_1min")
files.CAfr.preds3 <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.preds3)]
# Standardize data, required for PCA analysis
files.CAfr.preds3.scale <- scale(files.CAfr.preds3, center = TRUE, scale = TRUE)
# Calculate PCA for groups of correlated predictors.
pc.preds3 <- dudi.pca(files.CAfr.preds3.scale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pc.preds3$co # the column coordinates, which is the same as vector load for each PC
write.csv(pc.preds3$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_preds3_co_CAfr.csv")

# Define group 4 of collinear variables.
keep.preds4 <- c("bio14CAfr_1min", "bio17CAfr_1min", "bio2CAfr_1min", "bio19CAfr_1min")
files.CAfr.preds4 <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.preds4)]
# Standardize data, required for PCA analysis
files.CAfr.preds4.scale <- scale(files.CAfr.preds4, center = TRUE, scale = TRUE)
# Calculate PCA for groups of correlated predictors.
pc.preds4 <- dudi.pca(files.CAfr.preds4.scale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pc.preds4$co # the column coordinates, which is the same as vector load for each PC
write.csv(pc.preds4$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_preds4_co_CAfr.csv")

# Define group 5 of collinear variables.
keep.preds5 <- c("bio2CAfr_1min", "bio7CAfr_1min", "bio19CAfr_1min")
files.CAfr.preds5 <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.preds5)]
# Standardize data, required for PCA analysis
files.CAfr.preds5.scale <- scale(files.CAfr.preds5, center = TRUE, scale = TRUE)
# Calculate PCA for groups of correlated predictors.
pc.preds5 <- dudi.pca(files.CAfr.preds5.scale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pc.preds5$co # the column coordinates, which is the same as vector load for each PC
write.csv(pc.preds5$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_preds5_co_CAfr.csv")

# Define group 6 of collinear variables.
keep.preds6 <- c("DEM_maxmin_1min_CAfr", "DEM_sd_1min_CAfr")
files.CAfr.preds6 <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.preds6)]
# Standardize data, required for PCA analysis
files.CAfr.preds6.scale <- scale(files.CAfr.preds6, center = TRUE, scale = TRUE)
# Calculate PCA for groups of correlated predictors.
pc.preds6 <- dudi.pca(files.CAfr.preds6.scale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pc.preds6$co # the column coordinates, which is the same as vector load for each PC
write.csv(pc.preds6$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_preds6_co_CAfr.csv")

#######################################################################################
#########################  5. Preparing PCA axes as input variables for analysis  #####
#######################################################################################

# Based on the selection of environmental variables described above, a set of
# uncorrelated (Spearman rho < 0.7) is identified. On these predictors, we perform a
# PCA and export the first 2 PCA axes as .asc files.

# Select the uncorrelated environmental variables in 1 matrix.
keep.combined <- c("altCAfr_1min", "DEM_sd_1min_CAfr", "bio2CAfr_1min", "bio4CAfr_1min", "bio12CAfr_1min", "bio15CAfr_1min", "bio5CAfr_1min")
files.CAfr.combined.final.matrix <- files.CAfr.matrix[,(colnames(files.CAfr.matrix) %in% keep.combined)]

# Standardize data, required for PCA analysis.
files.CAfr.combined.final.matrix.scale <- scale(files.CAfr.combined.final.matrix, center = TRUE, scale = TRUE) 

# Perform the PCA using the dudi.pca() function.
pc.combined.final <- dudi.pca(files.CAfr.combined.final.matrix, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
# Get the vector loads for each PC.
pc.combined.final$co
write.csv(pc.combined.final$co, file = "C:/Mnspecies/Files/Multicollinearity/pc_co_CAfr.csv")
# Get the PC values for each raster cell.
PCexport.combined.final <- cbind(files.CAfr.df$y, files.CAfr.df$x, pc.combined.final$li)
colnames(PCexport.combined.final) <- c("y", "x", "PCA1", "PCA2", "PCA3", "PCA4")
PCexport.combined.final$PCA1 <- rescale(PCexport.combined.final$PCA1, to=c(0,1))
PCexport.combined.final$PCA2 <- rescale(PCexport.combined.final$PCA2, to=c(0,1))
summary(PCexport.combined.final)
PCexport.combined.final$PCA3 <- rescale(PCexport.combined.final$PCA3, to=c(0,1))
PCexport.combined.final$PCA4 <- rescale(PCexport.combined.final$PCA4, to=c(0,1))
summary(PCexport.combined.final)

# Below the adjusted function df2asc.
df2asc <- function (tdata, filenames = NULL, outdir = getwd(), gz = FALSE) 
{
  if (is.null(filenames)) {
    filenames = colnames(tdata)[3:length(tdata)]
  }
  else {
    if (length(filenames) != length(3:length(tdata))) 
      stop("variable names must be the same length as the files vector")
    filenames = as.character(filenames)
  }
  for (ii in 3:(length(tdata))) {
    lats = unique(tdata[, 1])
    lats = sort(lats)
    longs = unique(tdata[, 2])
    longs = sort(longs)
    cellsize = min(c(diff(lats), diff(longs)))
    # nc = ceiling((max(lats) - min(lats))/cellsize) + 1
    nc = ceiling((max(lats) - min(lats))/cellsize)
    # nr = ceiling((max(longs) - min(longs))/cellsize) + 1
    nr = ceiling((max(longs) - min(longs))/cellsize)
    out.asc = as.asc(matrix(NA, nr = nr, nc = nc), xll = min(longs), yll = min(lats), cellsize = cellsize)
    out.asc = put.data(tdata[, c(2:1, ii)], out.asc)
    write.asc(out.asc, paste(outdir, "/", filenames[ii - 2], sep = ""), gz = gz)
  }
}

# Export individual columns of the df as .asc files with the function df2asc(). # Note:
# the df has y, x coordinates (or lat,lon) and data columns in exacyly that exact order.
df2asc(PCexport.combined.final, outdir = "C:/Mnspecies/Files/Multicollinearity")

barplot(pc.combined.final$eig)
var1 <- pc.combined.final$eig[1]/sum(pc.combined.final$eig) # eigenvalue for PCA 1st axis
var1
var2 <- pc.combined.final$eig[2]/sum(pc.combined.final$eig) # eigenvalue for PCA 2nd axis
var2
var1+var2 # variance explained by first 2 axes
var3 <- pc.combined.final$eig[3]/sum(pc.combined.final$eig) # eigenvalue for PCA 3rd axis
var3
var4 <- pc.combined.final$eig[4]/sum(pc.combined.final$eig) # eigenvalue for PCA 4th axis
var4
var1 + var2 + var3 # variance explained by first 3 axes
var1 + var2 + var3 + var4 # variance explained by first 4 axes

#######################################################################################
#######################################################################################
#########################  END OF CODE  ###############################################
#######################################################################################
#######################################################################################
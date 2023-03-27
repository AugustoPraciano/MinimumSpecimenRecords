#######################################################################################
### Written by van Proosdij, A.S.J., Sosef, M.S.M., Wieringa, J.J. and Raes, N. 2016. 
### Adapted by Sampaio, A. C. P., Cavalcante, A. de M. B. 2021.
### Accurate Species Distribution Models: Minimum Required Number of Specimen Records
### in the Caatinga Biome
### Revista: / 
### Appendix 3: R script for Caatinga study area definition and analysis.
#######################################################################################

#######################################################################################
#########################  BRAZIL STUDY AREA SIMULATION  #############################
#######################################################################################

#######################################################################################
#######################################################################################
###  Written by André S.J. van Proosdij & Niels Raes, 2015
###  Adapted by Sampaio, A. C. P. (1) & Cavalcante, A. de M. B. (2), 2020
###  1 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
###  2 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
#######################################################################################
#######################################################################################

#######################################################################################
#########################  Index  #####################################################
#######################################################################################
# 1. Load functions and packages
# 2. Preparation of environmental variables
# 3. Settings of the analysis
# 4. Running the analysis
# 5. Null models
# 6. Summarize the results
# 7. Plotting the results
#######################################################################################

#######################################################################################
#########################  1. Load functions and packages  ############################
#######################################################################################

# Set the working directory. Create 5 folders: 'layers', 'outputs', 'outputsnm',
# 'presence', and 'projection'. Place the maxent.jar file in the library/dismo/java
# folder in C:/Program Files.
rm(list = ls(all = TRUE))
setwd("C:/Mnspecies/R/Realworld")
getwd()
dir.create("C:/Mnspecies/R/Realworld/layers")
dir.create("C:/Mnspecies/R/Realworld/outputs")
dir.create("C:/Mnspecies/R/Realworld/outputsnm")
dir.create("C:/Mnspecies/R/Realworld/presence")
dir.create("C:/Mnspecies/R/Realworld/projection")

library(raster) # for reading and writing rasters
library(stats) # for cor() function
library(mvtnorm) # for dmvnorm() function
library(SDMTools) # for write.asc(), asc2dataframe(), and pnt.in.poly() functions
library(dismo) # for evaluate() and maxent() functions
library(phyloclim) # for niche.overlap() function
library(plotrix) # for rescale() function
library(scales) # for rescale() function
library(ellipse) # for ellipse() function
library(ggplot2) # for ggplot() function
library(doBy)
library(grid)
removeTmpFiles(0)
source("C:/Mnspecies/R/sampaio_and_cavalcante_2020_A2_species_presence.R")
source("C:/Mnspecies/R/sampaio_and_cavalcante_2020_A4_nullmodel.R")

#######################################################################################
#########################  2. Preparation of environmental variables  #################
#######################################################################################

# All environmental layers need to be .asc files of exactly the same dimensions and
# should be in the base folder and in 'layers', 'presence', and 'projection' folders.
predictor.files <- list.files("C:/Mnspecies/R/Realworld/layers", pattern = '.asc', full.names = TRUE)
z <- stack(predictor.files) # Create a RasterStack
plot(z)


#######################################################################################
#########################  3. Settings of the analysis  ###############################
#######################################################################################

# The species prevalence is defined: the fraction of cells where the species is present.
prevalence.class <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Sample sizes with which the modelling will be done are defined.
sample_size <- c(3:20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 100)

# The number of repetitions is defined.
repetition <- c(1:100)

# The results of the analysis are placed in the vector 'joined_data'.
joined_data <- c()

# The predictor files are placed in a dataframe. For preparation of the PCA axes based on
# uncorrelated environmental variables, as well as preparation of the mask file used to
# select the ecological optimum of the simulated species see the script datapreparation.R
predictor.files <- list.files('C:/Mnspecies/R/Realworld/layers', pattern='.asc', full.names = TRUE)
predictor.names <- unlist(strsplit(basename(predictor.files), ".asc")) # Variable names
predictors.df <- asc2dataframe(predictor.files, varnames = predictor.names)

# The predictor files and the mask file are placed in a dataframe.
predictor.files.optimum <- list.files('C:/Mnspecies/R/Realworld', pattern='.asc', full.names = TRUE)
predictor.names.optimum <- unlist(strsplit(basename(predictor.files.optimum), ".asc")) # Variable names
predictors.df.optimum <- asc2dataframe(predictor.files.optimum, varnames = predictor.names.optimum)


#######################################################################################
#########################  4. Running the analysis  ###################################
#######################################################################################

for(w in repetition) {
  for(y in prevalence.class) {
    # A loop is used to converge the prevalence to predefined values. Prevalence depends
    # on the width of the ecological niche, which is set by the size of the standard
    # deviation of the species' bivariate normal function to 2 orthogonal environmental
    # variables. The loop starts with an initial, small value of sd ('start.sd'), which
    # is increased until the predefined prevalence is reached. The increase of start.sd
    # is relative to the difference between the realised and the desired prevalence.
    if (y == 0.1)
      prevalencethreshold <- 0.02
    if (y == 0.2)
      prevalencethreshold <- 0.01
    if (y == 0.3)
      prevalencethreshold <- 0.01
    if (y == 0.4)
      prevalencethreshold <- 0.005
    if (y == 0.5)
      prevalencethreshold <- 0.002
    
    #######################################################################################
    # Define the optimum values of the species respons to the environmental variables.
    #######################################################################################
    
    # A point is randomly selected from the Caatinga area by using the mask file.
    # The values of the environmental variables at this point are set as the mean values
    # of the bivariate normal response of the simulated species
    optimum.loc <- sample(nrow(predictors.df.optimum), 1) # Randomly sample 1 location
    optimum.values <- predictors.df.optimum[optimum.loc, ] # Extract predictor values
    
    # Plot to check for errors
    plot(z, 1, main = 'Random point defining species optimum'); points(optimum.values$x, optimum.values$y)
    
    keeps1 <- c("PCA1", "PCA2") # Keep the predictor values
    # 3D version: keeps1 <- c("PCA1", "PCA2", "PCA3")
    means <- as.vector(as.matrix(optimum.values[ ,(names(optimum.values) %in% keeps1)]))
    means <- as.numeric(as.character(means))
    
    # Start a loop to converge the prevalence towards the desired value of prevalence.
    start.sd <- 0.00001
    repeat {
      
      #######################################################################################
      # Define the habitat suitability and presence/absence of the simulated species by
      # running the species.presence() function.
      #######################################################################################
      
      # Sigma is the variance-covariance matrix of the bivariate normal distribution used to
      # define the habitat suitability. Sigma includes the standard deviation of each variable
      # (SD1 = SD2), covariance is 0 as the two variables are orthogonal.
      # Sigma = matrix(11, 12, 21, 22), with 11 = SD1, 12 = 21 = covarSD1*SD2, 22 = SD2.
      sigma <- matrix(
        data = c(start.sd, 0, 0, start.sd),
        nrow = 2,
        ncol = 2)
      # 3D version: sigma <- matrix(data = c(start.sd,0,0,0,start.sd,0,0,0,start.sd), nrow=3, ncol=3)
      
      # Run the species.presence() function on the environmental variables, means and sigma.
      true_presences <- species.presence(
        z     = z,
        means = means,
        sigma = sigma)
      
      # Retrieve the defined habitat suitability and presences from the species.presence output.
      suitability1 <- true_presences$suitability
      plot(suitability1)
      true_presences <- true_presences$presence
      
      # Convert the defined habitat suitability and true presence RasterLayers to vectors.
      suitability <- values(suitability1) # As vector, compare later with MaxEnt output
      pres <- na.omit(values(true_presences)) # Omit NA's
      
      # Calculate the realised prevalence of the species based on the SD of the predictors.
      prevalence <- length(pres[pres == 1]) / length(pres)
      
      # Check if the realised prevalence is aproximating the desired prevalence class value.
      h <- (y - prevalence)/y # Difference between realised prevalence and prevalence class
      if(h > prevalencethreshold | h < -prevalencethreshold)
        start.sd <- start.sd + 0.5*h*start.sd # Adjust start.sd and repeat the loop
      if(h >= -prevalencethreshold & h <= prevalencethreshold)
        break; # End the loop and continue with the next section
      
      # Save the defined presence/absence and habitat suitability for the simulated species.
      write.asc(
        x    = asc.from.raster(true_presences),
        file = "C:/Mnspecies/R/Realworld/presence/true_presences.asc");
      write.asc(
        x    = asc.from.raster(suitability1),
        file = "C:/Mnspecies/R/Realworld/presence/suitability.asc")
      
      # Objects to use in the rest of the analysis:
      #    prevalence: the realised prevalence value
      #    suitability: vector with defined habitat suitability
      #    true_presences: RasterLayer with defined presences and absences
    } # End of repeat function
    
    plot(true_presences) # Plot the given presences to check for errors
    
    for(x in sample_size) {
      
      # Read the predictor files, given suitability and given presences in one dataframe.
      # These are read all together from file to assure identical treatment of data.
      predictor.files2 <- list.files('C:/Mnspecies/R/Realworld/presence', pattern='.asc', full.names = TRUE)
      predictor.names2 <- unlist(strsplit(basename(predictor.files2), ".asc")) # Variable names
      predictors.df2 <- asc2dataframe(predictor.files2, varnames = predictor.names2)
      
      #######################################################################################
      # Sample the required number of records from the defined presences.
      #######################################################################################
      
      # Check if the number of given presences is larger than the required sample size.
      stopifnot(x < sum(predictors.df2$true_presences))
      
      # Add a column with ID numbers to identify rows.
      predictors.df2$ID <- c(1:nrow(predictors.df2))
      # Make a subset with given absences only.
      df2.absc <- predictors.df2[which(predictors.df2$true_presences == 0), ]
      # Make a subset with given presences only.
      df2.pres <- predictors.df2[which(predictors.df2$true_presences == 1), ]
      # Make a subset with sampled presences only. Sample the number of rows equal to the
      # sample size, habitat suitability ###TROCAR### is used as sampling probability.
      df2.presselect <- df2.pres[sample(nrow(df2.pres), x, prob = df2.pres$suitability), ]
      # Subtract the sampled presences from all given presences.
      df2.nonpres <- df2.pres[ !(df2.pres$ID %in% df2.presselect$ID), ]
      # Add the non-sampled given presences to the given absences.
      df2.pseudoabsc <- rbind(df2.nonpres, df2.absc)
      
      #######################################################################################
      # Run the function maxent() ('dismo' package, Hijmans & al., 2013)
      #######################################################################################
      
      # Create SWD files for presences and for background data as input for MaxEnt. Note:
      # By default, MaxEnt adds samples to background data.
      keeps2 <- c("PCA1", "PCA2") # Keep only the predictor variables
      pres <- df2.presselect[,(names(df2.presselect) %in% keeps2)] # Presence data
      absc <- df2.pseudoabsc[,(names(df2.pseudoabsc) %in% keeps2)] # Background data
      # Limit the number of background data to 10000.
      if (nrow(absc) > 10000) absc <- absc[sample(nrow(absc), 10000), ] else absc <- absc
      predictordata <- rbind(pres, absc) # rbind all presence and background data
      # Create a vector with 0/1 to identify presences and background.
      PAidentifier <- c(rep(1, nrow(pres)), rep(0, nrow(absc)))
      
      # Execute the maxent() function. By specifying the path, output files are stored. By
      # specifying the projection layers a prediction file "species_layers.asc" is generated.
      m <- maxent(predictordata,
                  PAidentifier,
                  path = "C:/Mnspecies/R/Realworld/outputs",
                  args = c("projectionlayers=layers", "redoifexists", "notooltips", "noautofeature", "linear", "quadratic", "nohinge", "noproduct", "nothreshold", "l2lqthreshold=1"))
      
      #######################################################################################
      # Get the MaxEnt AUC
      #######################################################################################
      
      # Read the MaxEnt training AUC value from the maxentResults file.
      maxentResults <- read.csv("C:/Mnspecies/R/Realworld/outputs/maxentResults.csv")
      MaxentAUC <- maxentResults$Training.AUC
      
      #######################################################################################
      # Calculate real AUC with the function evaluate() ('dismo' package, Hijmans & al., 2013)
      #######################################################################################
      
      # Real AUC is calculated using MaxEnt prediction values of given presences vs. given
      # absences. Data are in "predictedsuitability.v", "true_presences" identifies P/A's.
      realAUC.file <- 'C:/Mnspecies/R/Realworld/outputs/species_layers.asc'
      realAUC.names <- unlist(strsplit(basename(realAUC.file), ".asc")) # Variable names
      predictedsuitability.df <- asc2dataframe(realAUC.file, varnames = realAUC.names)
      
      # Concatenate objects with defined presences and absences and predicted suitability.
      allpoints <- cbind(predictors.df2, predictedsuitability.df)
      keeps3 <- c("true_presences", "species_layers")
      allpoints <- allpoints[,(names(allpoints) %in% keeps3)]
      pres.predictions <- allpoints[which(allpoints$true_presences == 1),2] # Select given P
      absc.predictions <- allpoints[which(allpoints$true_presences == 0),2] # Select given A
      # AUC functions can't handle very large data sets. For data sets > 10000, reduce the
      # number of presences and absences. This does not affect the result.
      if (length(absc.predictions) > 10000) {
        pres.predictions <- sample(pres.predictions, 0.1*length(pres.predictions))
        absc.predictions <- sample(absc.predictions, 0.1*length(absc.predictions))
      }
      AUCreal <- evaluate(p = pres.predictions, a = absc.predictions)@auc
      
      #######################################################################################
      # Calculate Spearman rank correlation using the function cor() ('stats' package,
      # R Core Team, 2014)
      #######################################################################################
      
      # Calculate the Spearman rank correlation of predicted vs. given habitat suitability.
      Spearman.files <- c("C:/Mnspecies/R/Realworld/outputs/species_layers.asc", "C:/Mnspecies/R/Realworld/presence/suitability.asc")
      Spearman.names <- unlist(strsplit(basename(Spearman.files), ".asc")) # Variable names
      predictors.Spearman <- asc2dataframe(Spearman.files, varnames = Spearman.names)
      predictedsuitability.v <- predictors.Spearman$species_layers
      givensuitability.v <- predictors.Spearman$suitability
      Spearman.cor <- cor(predictedsuitability.v, givensuitability.v, method = "spearman")
      
      #######################################################################################
      # Calculate niche overlap using the function niche.overlap() ('phyloclim' package, 
      # Heibl & Calenge, 2013)
      #######################################################################################
      
      #######################################################################################
      # Output
      #######################################################################################
      
      # The results are stored in a vector.
      a <- c(y, prevalence, x, w, start.sd, Spearman.cor, MaxentAUC, AUCreal, means[1], means[2])
      joined_data <- append(x = joined_data, values = a)
    }
  }
}


# Transform the vector 'joined_data' with output to a dataframe and name the columns.
joined_data <- matrix(
  data  = joined_data,
  byrow = TRUE,
  ncol  = 12)
joined_data <- data.frame(joined_data)

colnames(joined_data) <- c("prevalence", "true.prevalence", "sample.size", "repetition", "start.sd", "Spearman.cor", "Maxent.AUC", "AUCreal", "means1", "means2")


#######################################################################################
#########################  5 Null models  #############################################
#######################################################################################
# Run null models to test for significant deviance from chance (Raes & ter Steege, 2007).

# Get the unique sample sizes.
joined_data2 <- joined_data
records <- unique(joined_data2$sample.size)
records

# Prepare the environmental data.
x <- predictors.df # Predictor files as created above
drops2 <- c("y", "x", "true_presences", "suitability")
x <- x[,!(names(x) %in% drops2)]

# Run the loop for all unique sample sizes.
nullAUC <- list() # Create an empty list to store the results
for (n in records) {
  a <- nullModelreal(x, n, rep = 9)
  nullAUC[[n]] <- a # Store the null AUC values for each sample size
}

# Replace NULL with NA for sample sizes that are not present.
for (i in 1:length(nullAUC)) {
  if (is.null(nullAUC[[i]]) == TRUE) {nullAUC[[i]] <- NA}
}

# Turn nullAUC into a data.frame with all AUC values for a sample size in one column.
nullAUC <- do.call(cbind, nullAUC)
nullAUCdataframe <- as.data.frame(nullAUC)

# For each model, get the rank number of the real AUC compared to null model AUC's.
rankAUCreal <- list() # Create an empty vector to store everything
joined_data2 <- as.matrix(joined_data)
for (i in 1:nrow(joined_data2)) { # For each unique repetition...
  j <- as.numeric(as.character(joined_data2[i,3])) # ...get the sample size
  k <- nullAUCdataframe[,j] # ...get the null model AUC values of that sample size
  l <- joined_data2[i,8] # ...get the AUCreal value of that sample size
  m <- c(as.numeric(as.character(l)), as.numeric(as.character(k))) # ...place all in 1 list
  n <- rank(m) # ...rank them
  o <- n[1] # ...get the rank of the AUCreal
  rankAUCreal[i] <- o # ...and store the rank of AUCreal
}
rankAUCrealunlist <- unlist(rankAUCreal) # Unlist rankAUCreal before cbinding
joined_data3 <- cbind(joined_data2, rankAUCrealunlist) # Attach the column with rank numbers
colnames(joined_data3)[13] <- "rankAUCreal" # Name the column

# For each model, get the rank number of the AUCMaxEnt compared to null model AUC's.
rankAUCMaxent <- list() # Create an empty vector to store everything
for (i in 1:nrow(joined_data3)) { # For each unique repetition...
  j <- as.numeric(as.character(joined_data3[i,3])) # ...get the sample size
  k <- nullAUCdataframe[,j] # ...get the null model AUC values of that sample size
  l <- joined_data3[i,7] # ...get the AUCMaxEnt value of that sample size
  m <- c(as.numeric(as.character(l)), as.numeric(as.character(k))) # ...place all in 1 list
  n <- rank(m) # ...rank them
  o <- n[1] # ...get the rank of the AUCMaxEnt
  rankAUCMaxent[i] <- o # ...and store the rank of AUCMaxEnt
}
rankAUCMaxentunlist <- unlist(rankAUCMaxent) # Unlist rankAUCMaxent before cbinding
joined_data4 <- cbind(joined_data3, rankAUCMaxentunlist) # Attach the column with rank numbers
colnames(joined_data4)[14] <- "rankAUCMaxent" # Name the column

# For each model calculate the difference between the rank of real AUC and MaxEnt AUC.
diff.rank.AUC.real.Maxent <- list() # Create an empty vector to store everything
for (i in 1:nrow(joined_data4)) { # For each unique repetition...
  p <- as.numeric(as.character(joined_data4[i,13])) # ...load the rank of the real AUC
  q <- as.numeric(as.character(joined_data4[i,14])) # ...load the rank of the MAxEnt AUC
  diff.rank.AUC.real.Maxent[i] <- p-q # ... and calculate the difference
}
diff.rank.AUC.real.Maxent.unlist <- unlist(diff.rank.AUC.real.Maxent) # Unlist before cbinding it
joined_data5 <- cbind(joined_data4, diff.rank.AUC.real.Maxent.unlist) # Attach the column with the difference in rank numbers between the AUCreal and AUCMaxent
colnames(joined_data5)[15] <- "diffrankAUCrealMaxent"


#######################################################################################
#########################  6 Summarize the results  ###################################
#######################################################################################

# For each prevalence class and each sample size, get the lower and upper limit of the
# upper 95% range of the values. This effectively excludes the 5% worst performing models.

# Remove all rows which contain NA's.
joined_data6 <- as.data.frame(joined_data5)
joined_data6 <- na.omit(joined_data6)

# Spearman rank correlation
summarySpearmanreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- joined_data6[(joined_data6$prevalence == i), ]
  for (j in sample_size) { # For each sample size
    abc2 <- abc[(abc$sample.size == j), ]
    datasort <- sort(abc2$Spearman.cor) # Sort the values increasing
    lowerlimit <- datasort[6] # Get the lower limit of the upper 95% range of values
    upperlimit <- max(datasort) # Get the upper limit
    abc3 <- cbind(i, j, lowerlimit, upperlimit) # Prevalence, sample size, lower and upper limit
    summarySpearmanreal <- rbind(summarySpearmanreal, abc3)
  }
}
colnames(summarySpearmanreal) <- c("prevalence", "sample.size", "Spearmanll", "Spearmanul")
write.csv(x = summarySpearmanreal, file = "C:/Mnspecies/R/Realworld/summarySpearmanreal.csv")

# MaxEnt AUC values
summaryMaxEntAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- joined_data6[(joined_data6$prevalence == i), ]
  for (j in sample_size) { # For each sample size
    abc2 <- abc[(abc$sample.size == j), ]
    datasort <- sort(abc2$Maxent.AUC) # Sort the values increasing
    lowerlimit <- datasort[6] # Get the lower limit of the upper 95% range of values
    upperlimit <- max(datasort) # Get the upper limit
    abc3 <- cbind(i, j, lowerlimit, upperlimit) # Prevalence, sample size, lower and upper limit
    summaryMaxEntAUCreal <- rbind(summaryMaxEntAUCreal, abc3)
  }
}
colnames(summaryMaxEntAUCreal) <- c("prevalence", "sample.size", "MaxEntAUCll", "MaxEntAUCul")
write.csv(x = summaryMaxEntAUCreal, file = "C:/Mnspecies/R/Realworld/summaryMaxEntAUCreal.csv")

# Real AUC values
summaryrealAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- joined_data6[(joined_data6$prevalence == i), ]
  for (j in sample_size) { # For each sample size
    abc2 <- abc[(abc$sample.size == j), ]
    datasort <- sort(abc2$AUCreal) # Sort the values increasing
    lowerlimit <- datasort[6] # Get the lower limit of the upper 95% range of values
    upperlimit <- max(datasort) # Get the upper limit
    abc3 <- cbind(i, j, lowerlimit, upperlimit) # Prevalence, sample size, lower and upper limit
    summaryrealAUCreal <- rbind(summaryrealAUCreal, abc3)
  }
}
colnames(summaryrealAUCreal) <- c("prevalence", "sample.size", "realAUCll", "realAUCul")
write.csv(x = summaryrealAUCreal, file = "C:/Mnspecies/R/Realworld/summaryrealAUCreal.csv")

# MaxEnt AUC rank values
summaryrankMaxEntAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- joined_data6[(joined_data6$prevalence == i), ]
  for (j in sample_size) { # For each sample size
    abc2 <- abc[(abc$sample.size == j), ]
    datasort <- sort(abc2$rankAUCMaxent) # Sort the values increasing
    lowerlimit <- datasort[6] # Get the lower limit of the upper 95% range of values
    upperlimit <- max(datasort) # Get the upper limit
    abc3 <- cbind(i, j, lowerlimit, upperlimit) # Prevalence, sample size, lower and upper limit
    summaryrankMaxEntAUCreal <- rbind(summaryrankMaxEntAUCreal, abc3)
  }
}
colnames(summaryrankMaxEntAUCreal) <- c("prevalence", "sample.size", "rankMaxEntAUCll", "rankMaxEntAUCul")
write.csv(x = summaryrankMaxEntAUCreal, file = "C:/Mnspecies/R/Realworld/summaryrankMaxEntAUCreal.csv")

# Real AUC rank values
summaryrankrealAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- joined_data6[(joined_data6$prevalence == i), ]
  for (j in sample_size) { # For each sample size
    abc2 <- abc[(abc$sample.size == j), ]
    datasort <- sort(abc2$rankAUCreal) # Sort the values increasing
    lowerlimit <- datasort[6] # Get the lower limit of the upper 95% range of values
    upperlimit <- max(datasort) # Get the upper limit
    abc3 <- cbind(i, j, lowerlimit, upperlimit) # Prevalence, sample size, lower and upper limit
    summaryrankrealAUCreal <- rbind(summaryrankrealAUCreal, abc3)
  }
}
colnames(summaryrankrealAUCreal) <- c("prevalence", "sample.size", "rankrealAUCll", "rankrealAUCul")
write.csv(x = summaryrankrealAUCreal, file = "C:/Mnspecies/R/Realworld/summaryrankrealAUCreal.csv")

#######################################################################################
# Fitting curves on the data.
#######################################################################################
# To mask out small stochastic effects, per prevalence class, smooth the lower and upper
# limits of the upper 95% range of the model performance values using the loess() function
# ('stats' package, R Core Team, 2014).

# Spearman rank values.
summarySpearmanreal <- as.data.frame(summarySpearmanreal)
summarySpearmanreal2 <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summarySpearmanreal[(summarySpearmanreal$prevalence == i), ]
  Spearmanll.loess <- loess(Spearmanll ~ sample.size, data = abc) # Fit a curve
  Spearmanul.loess <- loess(Spearmanul ~ sample.size, data = abc) # Fit a curve
  Spearmanllsmooth <- predict(Spearmanll.loess) # Get the smoothed values
  Spearmanulsmooth <- predict(Spearmanul.loess) # Get the smoothed values
  abc2 <- cbind(abc$prevalence, abc$sample.size, Spearmanllsmooth, Spearmanulsmooth)
  summarySpearmanreal2 <- rbind(summarySpearmanreal2, abc2)
}
colnames(summarySpearmanreal2) <- c("prevalence", "sample.size", "Spearmanllsmooth", "Spearmanulsmooth")
summarySpearmanreal2 <- as.data.frame(summarySpearmanreal2)
summarySpearmanreal2$Spearmanllsmooth[summarySpearmanreal2$Spearmanllsmooth > 1] <- 1 # Truncate all Spearmanllsmooth values > 1 to 1
summarySpearmanreal2$Spearmanulsmooth[summarySpearmanreal2$Spearmanulsmooth > 1] <- 1 # Truncate all Spearmanulsmooth values > 1 to 1

# MaxEnt AUC values.
summaryMaxEntAUCreal <- as.data.frame(summaryMaxEntAUCreal)
summaryMaxEntAUCreal2 <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryMaxEntAUCreal[(summaryMaxEntAUCreal$prevalence == i), ]
  MaxEntAUCll.loess <- loess(MaxEntAUCll ~ sample.size, data = abc) # Fit a curve
  MaxEntAUCul.loess <- loess(MaxEntAUCul ~ sample.size, data = abc) # Fit a curve
  MaxEntAUCllsmooth <- predict(MaxEntAUCll.loess) # Get the smoothed values
  MaxEntAUCulsmooth <- predict(MaxEntAUCul.loess) # Get the smoothed values
  abc2 <- cbind(abc$prevalence, abc$sample.size, MaxEntAUCllsmooth, MaxEntAUCulsmooth)
  summaryMaxEntAUCreal2 <- rbind(summaryMaxEntAUCreal2, abc2)
}
colnames(summaryMaxEntAUCreal2) <- c("prevalence", "sample.size", "MaxEntAUCllsmooth", "MaxEntAUCulsmooth")
summaryMaxEntAUCreal2 <- as.data.frame(summaryMaxEntAUCreal2)
summaryMaxEntAUCreal2$MaxEntAUCllsmooth[summaryMaxEntAUCreal2$MaxEntAUCllsmooth > 1] <- 1 # Truncate all MaxEntAUCllsmooth values > 1 to 1
summaryMaxEntAUCreal2$MaxEntAUCulsmooth[summaryMaxEntAUCreal2$MaxEntAUCulsmooth > 1] <- 1 # Truncate all MaxEntAUCulsmooth values > 1 to 1

# Real AUC values.
summaryrealAUCreal <- as.data.frame(summaryrealAUCreal)
summaryrealAUCreal2 <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryrealAUCreal[(summaryrealAUCreal$prevalence == i), ]
  realAUCll.loess <- loess(realAUCll ~ sample.size, data = abc) # Fit a curve
  realAUCul.loess <- loess(realAUCul ~ sample.size, data = abc) # Fit a curve
  realAUCllsmooth <- predict(realAUCll.loess) # Get the smoothed values
  realAUCulsmooth <- predict(realAUCul.loess) # Get the smoothed values
  abc2 <- cbind(abc$prevalence, abc$sample.size, realAUCllsmooth, realAUCulsmooth)
  summaryrealAUCreal2 <- rbind(summaryrealAUCreal2, abc2)
}
colnames(summaryrealAUCreal2) <- c("prevalence", "sample.size", "realAUCllsmooth", "realAUCulsmooth")
summaryrealAUCreal2 <- as.data.frame(summaryrealAUCreal2)
summaryrealAUCreal2$realAUCllsmooth[summaryrealAUCreal2$realAUCllsmooth > 1] <- 1 # Truncate all realAUCllsmooth values > 1 to 1
summaryrealAUCreal2$realAUCulsmooth[summaryrealAUCreal2$realAUCulsmooth > 1] <- 1 # Truncate all realAUCulsmooth values > 1 to 1

# MaxEnt AUC rank values.
summaryrankMaxEntAUCreal <- as.data.frame(summaryrankMaxEntAUCreal)
summaryrankMaxEntAUCreal2 <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryrankMaxEntAUCreal[(summaryrankMaxEntAUCreal$prevalence == i), ]
  rankMaxEntAUCll.loess <- loess(rankMaxEntAUCll ~ sample.size, data = abc) # Fit a curve
  rankMaxEntAUCul.loess <- loess(rankMaxEntAUCul ~ sample.size, data = abc) # Fit a curve
  rankMaxEntAUCllsmooth <- predict(rankMaxEntAUCll.loess) # Get the smoothed values
  rankMaxEntAUCulsmooth <- predict(rankMaxEntAUCul.loess) # Get the smoothed values
  abc2 <- cbind(abc$prevalence, abc$sample.size, rankMaxEntAUCllsmooth, rankMaxEntAUCulsmooth)
  summaryrankMaxEntAUCreal2 <- rbind(summaryrankMaxEntAUCreal2, abc2)
}
colnames(summaryrankMaxEntAUCreal2) <- c("prevalence", "sample.size", "rankMaxEntAUCllsmooth", "rankMaxEntAUCulsmooth")
summaryrankMaxEntAUCreal2 <- as.data.frame(summaryrankMaxEntAUCreal2)
summaryrankMaxEntAUCreal2$rankMaxEntAUCllsmooth[summaryrankMaxEntAUCreal2$rankMaxEntAUCllsmooth > 100] <- 100 # Truncate all rankMaxEntAUCllsmooth values > 100 to 100
summaryrankMaxEntAUCreal2$rankMaxEntAUCulsmooth[summaryrankMaxEntAUCreal2$rankMaxEntAUCulsmooth > 100] <- 100 # Truncate all rankMaxEntAUCulsmooth values > 100 to 100

# Real AUC rank values.
summaryrankrealAUCreal <- as.data.frame(summaryrankrealAUCreal)
summaryrankrealAUCreal2 <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryrankrealAUCreal[(summaryrankrealAUCreal$prevalence == i), ]
  rankrealAUCll.loess <- loess(rankrealAUCll ~ sample.size, data = abc) # Fit a curve
  rankrealAUCul.loess <- loess(rankrealAUCul ~ sample.size, data = abc) # Fit a curve
  rankrealAUCllsmooth <- predict(rankrealAUCll.loess) # Get the smoothed values
  rankrealAUCulsmooth <- predict(rankrealAUCul.loess) # Get the smoothed values
  abc2 <- cbind(abc$prevalence, abc$sample.size, rankrealAUCllsmooth, rankrealAUCulsmooth)
  summaryrankrealAUCreal2 <- rbind(summaryrankrealAUCreal2, abc2)
}
colnames(summaryrankrealAUCreal2) <- c("prevalence", "sample.size", "rankrealAUCllsmooth", "rankrealAUCulsmooth")
summaryrankrealAUCreal2 <- as.data.frame(summaryrankrealAUCreal2)
summaryrankrealAUCreal2$rankrealAUCllsmooth[summaryrankrealAUCreal2$rankrealAUCllsmooth > 100] <- 100 # Truncate all rankrealAUCllsmooth values > 100 to 100
summaryrankrealAUCreal2$rankrealAUCulsmooth[summaryrankrealAUCreal2$rankrealAUCulsmooth > 100] <- 100 # Truncate all rankrealAUCulsmooth values > 100 to 100

#######################################################################################
# Export the summarized results.
#######################################################################################

write.csv(x = summarySpearmanreal2, file = "C:/Mnspecies/R/Realworld/summarySpearmanreal2.csv")
write.csv(x = summaryMaxEntAUCreal2, file = "C:/Mnspecies/R/Realworld/summaryMaxEntAUCreal2.csv")
write.csv(x = summaryrealAUCreal2, file = "C:/Mnspecies/R/Realworld/summaryrealAUCreal2.csv")
write.csv(x = summaryrankMaxEntAUCreal2, file = "C:/Mnspecies/R/Realworld/summaryrankMaxEntAUCreal2.csv")
write.csv(x = summaryrankrealAUCreal2, file = "C:/Mnspecies/R/Realworld/summaryrankrealAUCreal2.csv")


#######################################################################################
#########################  7. Read the minimum required sample sizes  #################
#######################################################################################

# Identify the minimum sample size for which the lower limit of the upper 95% range of
# the model performance values exceeds the defined critical value.

# Lower limit of the upper 95% range of Spearman rank correlation values > 0.9.
summarySpearmanreal2 <- as.data.frame(summarySpearmanreal2)
minimumrecordsSpearmanreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summarySpearmanreal2[(summarySpearmanreal2$prevalence == i), ]
  dataselect <- abc[(abc$Spearmanllsmooth > 0.9), ]
  criticalsamplesize <- min(dataselect$sample.size)
  abc2 <- c(i, criticalsamplesize)
  minimumrecordsSpearmanreal <- rbind(minimumrecordsSpearmanreal, abc2)
}
colnames(minimumrecordsSpearmanreal) <- c("prevalence", "minimumrecords")
minimumrecordsSpearmanreal

# Lower limit of the upper 95% range of real AUC values > 0.9.
summaryrealAUCreal2 <- as.data.frame(summaryrealAUCreal2)
minimumrecordsrealAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryrealAUCreal2[(summaryrealAUCreal2$prevalence == i), ]
  dataselect <- abc[(abc$realAUCllsmooth > 0.9), ]
  criticalsamplesize <- min(dataselect$sample.size)
  abc2 <- c(i, criticalsamplesize)
  minimumrecordsrealAUCreal <- rbind(minimumrecordsrealAUCreal, abc2)
}
colnames(minimumrecordsrealAUCreal) <- c("prevalence", "minimumrecords")
minimumrecordsrealAUCreal

# Lower limit of the upper 95% range of MaxEnt AUC rank values > 95.
summaryrankMaxEntAUCreal2 <- as.data.frame(summaryrankMaxEntAUCreal2)
minimumrecordsrankMaxEntAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryrankMaxEntAUCreal2[(summaryrankMaxEntAUCreal2$prevalence == i), ]
  dataselect <- abc[(abc$rankMaxEntAUCllsmooth > 95), ]
  criticalsamplesize <- min(dataselect$sample.size)
  abc2 <- c(i, criticalsamplesize)
  minimumrecordsrankMaxEntAUCreal <- rbind(minimumrecordsrankMaxEntAUCreal, abc2)
}
colnames(minimumrecordsrankMaxEntAUCreal) <- c("prevalence", "minimumrecords")
minimumrecordsrankMaxEntAUCreal

# Lower limit of the upper 95% range of real AUC rank values > 95.
summaryrankrealAUCreal2 <- as.data.frame(summaryrankrealAUCreal2)
minimumrecordsrankrealAUCreal <- c()
for (i in prevalence.class) { # For each prevalence class
  abc <- summaryrankrealAUCreal2[(summaryrankrealAUCreal2$prevalence == i), ]
  dataselect <- abc[(abc$rankrealAUCllsmooth > 95), ]
  criticalsamplesize <- min(dataselect$sample.size)
  abc2 <- c(i, criticalsamplesize)
  minimumrecordsrankrealAUCreal <- rbind(minimumrecordsrankrealAUCreal, abc2)
}
colnames(minimumrecordsrankrealAUCreal) <- c("prevalence", "minimumrecords")
minimumrecordsrankrealAUCreal

#######################################################################################
#######################################################################################
#########################  END OF CODE  ###############################################
#######################################################################################
#######################################################################################
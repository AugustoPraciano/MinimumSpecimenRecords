#######################################################################################
### written by van Proosdij, A.S.J., Sosef, M.S.M., Wieringa, J.J. and Raes, N. 2016. Adapted by
### Sampaio, A. C. P., Cavalcante, A. de M. B. 2021.
### Accurate Species Distribution Models: Minimum Required Number of Specimen Records
### in the Caatinga Biome
### Revista: / 
### Appendix 2: R script for simulated species definition for the Caatinga study area.
#######################################################################################

#######################################################################################
#########################  SPECIES.PRESENCE FUNCTION  #################################
#######################################################################################

#######################################################################################
#######################################################################################
###  Written by André S.J. van Proosdij, & Niels Raes, 2015 Adapted by
###  Sampaio, A. C. P. (1), Cavalcante, A. de M. B. (2) 2020.
###  1 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
###  2 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
#######################################################################################
#######################################################################################

# DESCRIPTION
# The 'species.presence' function creates a habitat suitability and a presence/absence
# map. Habitat suitability is defined as a multivariate normal function of the
# environmental variables using the function dmvnorm() ('mvtnorm' package, Genz & al.,
# 2014).
# To generate presence/absence map, for each point is checked whether it is inside the
# parameter space defined by the species multivariate respons to the variables. For
# this, first the variables are converted into independent standard normal variables
# using Choleski Decomposition implemented in the function chol(). Then, for each point
# the chi-squared distance to the origin is calculated using the function qchisq()
# ('stats' package, R Core Team, 2014), with the number of dimensions as degrees of
# freedom. Points which distance does not exceed the threshold value are defined as
# presence, other points are defined as absences. Here, the function is applied to
# analysis using 2 variables, but it can be applied to more dimensions too.

# ARGUMENTS
# z       - The environmental variables (predictors) that are used to define the habitat
#           suitability of each raster cell. This argument should be a RasterStack with
#           RasterLayers for each environmental variable.
# means   - The mean values for the variables of the bivariate normal distribution. This
#           argument should be a vector.
# sigma   - The variance-covariance matrix of the bivariate normal distribution. This
#           argument should be a matrix.

species.presence <- function(z, means, sigma) {
  #######################################################################################
  #########################  Calculate the density function #############################
  #######################################################################################
  
  # Transform the RasterStack with predictors ('z') to a matrix, required by dmvnorm().
  z1 <- as.matrix(as.data.frame(z))
  # Calculate the habitat suitabiliy as function of 'z', 'means' and 'sigma' using the
  # function dmvnorm() ('mvtnorm' package).
  suitability <- dmvnorm(x = z1, mean = means, sigma = sigma)
  
  #######################################################################################
  #########################  Calculate the presence and absence of the species  #########
  #######################################################################################
  
  # Define the number of dimensions and threshold level defining presences and absences.
  dim <- ncol(z1)
  level <- 0.68
  
  # Convert the variabes into independent standard normal variables using Choleski
  # Decomposition implemented in the function chol().
  # First, set the mean values to '0' and get the distance of each point to the means.
  compmeans <- matrix(means, nrow = nrow(z1), ncol = dim, byrow = TRUE)
  pc <- z1 - compmeans
  # Construct the covariance unit.
  M <- t(chol(sigma))
  ip <- t(solve(M) %*% t(pc))

  # Determine B such that (rowSums(ip^2) < B) = level. Use the chi-squared distribution
  # implemented in the function qchisq() ('stats' package, R Core Team, 2014), with the
  # number of dimensions as degrees of freedom.
  B <- qchisq(level, df = dim)
  # Check for each point if it is located inside or outside the object.
  inellipse <- rowSums(ip^2) < B
  presence <- cbind(z1, inellipse)
  presence <- presence[ ,"inellipse"]
  
  # When applied to 2 variabbles (dim = 2), view the result.
  if(dim == 2)
  {
    plot(z1[ ,1], z1[ ,2])
    ellipse2 <- ellipse(sigma, centre = means, level = level)
    lines(ellipse2, col = 'red')
  }
  
  #######################################################################################
  #########################  Output  ####################################################
  #######################################################################################
  
  # The output consists of a list containing a RasterLayer for the habitat suitability
  # and a RasterLayer containing presences and absences represented by 1's and 0's.
  r.presence <- r.suitability <- raster(z) # Create two rasters
  values(r.suitability) <- suitability # Fill the raster with habitat suitability values
  values(r.presence) <- presence # Fill the raster with presence/absence values
  out.list <- list(r.suitability, r.presence)
  names(out.list) <- c("suitability", "presence")
  return(out.list)
}

#######################################################################################
#########################  END OF CODE  ###############################################
#######################################################################################
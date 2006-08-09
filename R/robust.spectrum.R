### robust.spectrum.R  (2005-05-24)
###
###    Robust estimate of spectrum of time series
###
### Copyright 2005 Miika Ahdesmaki
###
###
###
### This file is part of the `GeneCycle' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA




######################################################################
######################################################################
# robust.spectrum.R
######################################################################
######################################################################

# public functions


#################################################################
############# The spectra of multiple time series ###############
#################################################################
robust.spectrum <- function(x) 
{
  # x should be a matrix consisting of the time series as column
  # vectors.
  if (is.longitudinal(x))
  {
    if (has.repeated.measurements(x) )
      stop("Multiple measurements per time point are not supported")
    
    if (!is.equally.spaced(x))
      warning("Function assumes equally spaced time points")
  }
  else
  {
    x <- as.matrix(x)
  }

  noTS <- ncol(x) #number of time series
  y <- matrix(NA,nrow(x),ncol(x))
  for (h in 1:noTS)
  {
    temp <- x[,h]
    y[,h] <- robust.spectrum.single(temp)
  }
  
  return(y)
}

######################################################################

# private functions


#################################################################
# This function implements the robust, rank-based spectral#######
# estimator introduced in Pearson et al. 2003####################
#################################################################
robust.spectrum.single <- function(x) 
{
  
  if( is.constant.single(x) ) warning("Constant time series!")
  
  ##############################################
  # Some adjustable parameters
  ##############################################
  # Length of the zero-padded one-sided "rho"(=Rsm)
  zp <- 2*length(x)
  
  # Length of the original sequence
  n <- length(x)
  
  # Let us define the maximum lag for the correlation coefficient:
  maxM <- n-2
  
  
  ##############################################
  # Correlation coefficient
  ##############################################
  # Reserve space
  Rsm <- matrix(NA, nrow = maxM+1, ncol = 1)
    
  
  missing <- is.na(x)
  nonMissing <- !missing
  mi <- which(missing)     # indices for missing values
  nmi <- which(nonMissing) # indices for nonmissing values
  
  # Mean removal
  x[nmi] <- x[nmi] - mean(x[nmi])
  
  # Run through all the lags
  for (lags in 0:maxM)
  {
    # Modified Spearman's method
    indexes <- 1:(n-lags)	# Initial indices
    ends <- length(nonMissing)
  
    # Values in both the original and shifted vectors must be present:
    temp <- (nonMissing[1:(ends-lags)] + nonMissing[(lags+1):ends]) >= 2
    indexPresent <- which(temp)
    indexes <- indexes[indexPresent]	# The indices that are present in
  							# both sequences
    Rsm[lags+1] <- ifelse(
                     length(indexes)<=1 , 
		     0 ,
		     spearman(x[indexes],x[indexes+lags],n))
  }
  
  # Zero-padding
  Rsm[(length(Rsm)+1):zp] <- 0
  fftemp <- fft(Rsm)
  
  # The following implementation is as in (Ahdesmäki, Lähdesmäki et al., 2005)
  Ssm <- abs( 2*Re(fftemp) - Rsm[1] )
  Ssm <- Ssm[1:floor(length(Ssm)/2)]
  
  # Return the spectral content, frequencies [0,pi)
  return(Ssm)
}


#################################################################
########Inner function:Spearman's correlation coefficient########
#################################################################


spearman <- function(x, y, N, version=c("builtin", "miika") )
{
  #cat (paste("DEBUG: ", length(x), length(y), N, "\n") )
  
  version <- match.arg(version)
  
  if (version == "builtin")
  {
      rho <- cor(x, y, method="spearman" ) * length(x)/N
  }
  
  if (version == "miika")
  {
      Km <- length(x)
      sx <- sort(x); ix <- order(x);
      sy <- sort(y); iy <- order(y);

      rx <- matrix(0, 1, Km)
      ry <- rx
      rx[ix] <- (1:Km)
      ry[iy] <- (1:Km)

      rho <- 12/(N*(Km^2-1)) * (rx-(Km+1)/2) %*% (t(ry) - (Km+1)/2)
  }
  
  return(rho)
}


#############################################
# check that the two versions are the same
#
# x1 <- rnorm(10)
# x2 <- rnorm(10)
#
# spearman(x1, x2, 111)
# spearman(x1, x2, 111, version="miika")
#############################################


######################################################################

\name{robust.g.test}
\alias{robust.g.test}
\alias{g.statistic}
\alias{robust.spectrum}


\title{Robust g Test for Multiple (Genetic) Time Series}
\description{
  \code{robust.g.test} calculates the p-value(s) for a robust
   nonparametric version of Fisher's g-test (1929).  Details
  of this approach are described in Ahdesmaki et al. (2005), along with
  an extensive discussion of its application to gene expression data.
  
  \code{g.statistic} computes the test statistic given a discrete
   time series spectrum.
  
  \code{robust.spectrum} computes a robust rank-based estimate
  of the periodogram/correlogram - see  Ahdesmaki et al. (2005)
  for details.
}
\usage{
robust.g.test(y, index, perm = FALSE, x, noOfPermutations = 5000)
g.statistic(y, index)
robust.spectrum(x)
}
\arguments{
  \item{y}{the matrix consisting of the spectral estimates
            as column vectors}
  \item{x}{a matrix consisting of the time series as column
           vectors.  In \code{robust.g.test} only needed if
           permutation  tests are used}
  \item{index}{an index to the spectral
               estimates that is to be used in the testing for periodicity.
               If \code{index} is missing, the maximum component of the
	       spectral estimate is used in testing (regardless of the
	       frequency of this maximum)}
  \item{perm}{if \code{perm} is \code{FALSE}, 
              a simulated distribution for the g-statistic is used. 
              If per \code{perm} is \code{TRUE}, permutation  tests are 
	      used to find  the distribution of the g-statistic for each 
	      time series separately. }
  
  \item{noOfPermutations}{number of permutations that 
                          are used for each time series 
			  (default = 5000)}
}
\details{
  Application of \code{robust.g.test} can be very computer intensive, especially
  the production of the distribution of the test statistics may take a lot
  of time. Therefore, this distribution (dependening on the length of
  the time series) is stored in an external file to avoid recomputation
  (see example below).
  
  For the general idea behind the Fisher's g test also see  
  \code{\link{fisher.g.test}} which implements an analytic approach for g-testing.
  This is faster but not robust and also assumes Gaussian noise.
  
}
\value{
 \code{robust.g.test} returns a list of p-values, and \code{g.statistic} the associated
 test statistics.  \code{robust.spectrum} returns a matrix where the column vectors correspond
 to the spectra corresponding to each time series.
}
\author{
  Miika Ahdesmaki (\email{miika.ahdesmaki@tut.fi}).
}
\references{
  Fisher, R.A. (1929). Tests of significance in harmonic analysis. \emph{Proc. Roy. Soc. A}, \bold{125}, 54--59.


  Ahdesmaki, M., Lahdesmaki, H., Peason, R., Huttunen, H., and Yli-Harja O.
  (2005).  BMC Bioinformatics 6:117.


}
\seealso{\code{\link[fdrtool]{fdrtool}}, \code{\link{fisher.g.test}}.}

\examples{
\dontrun{

# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# how many samples and and how many genes?
dim(caulobacter)


# robust, rank-based spectral estimator applied to first 5 genes
spe5 <- robust.spectrum(caulobacter[,1:5])

# g statistics computed from the spectrum
g.statistic(spe5)

# robust p-values
pval <- robust.g.test(spe5)  # generates a file with the name "g_pop_length_11.txt"
pval <- robust.g.test(spe5)  # second call: much faster..

pval

# delete the external file 
unlink("g_pop_length_11.txt") 

}
}
\keyword{htest}

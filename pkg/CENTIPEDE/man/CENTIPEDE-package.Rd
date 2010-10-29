\name{CENTIPEDE-package}
\alias{CENTIPEDE-package}
\alias{CENTIPEDE}
\alias{Centipede}
\alias{centipede}
\docType{package}

\title{
CENTIPEDE a probabilistic model for learning DNaseI footprints and
predicting TF binding sites
}

\description{
  Centipede fits a bayesian hierarchical mixture model to
  learn TF-specific distribution of experimental data on a particular 
  cell-type for a set of candidate binding sites described by a motif. 
  More documentation is under preparation and will be make available
  soon at http://centipede.uchicago.edu
}

\details{
\tabular{ll}{
Package: \tab CENTIPEDE\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2010-10-27\cr
License: \tab GPL \cr
URL: \tab http://centipede.uchicago.edu \cr
LazyLoad: \tab yes\cr
}
}

\author{
  Roger Pique-Regi <rpique@gmail.com> and Jacob F Degner <jdegner@uchicago.edu>
  Maintainer: Roger Pique-Regi <rpique@uchicago.edu>
}


\references{
  Roger Pique-Regi, Jacob F. Degner, Athma A. Pai, Daniel J. Gaffney,
  Yoav Gilad, Jonathan K. Pritchard. "Accurate inference of
  transcription factor binding from DNA sequence and chromatin
  accessibility data", Genome Research, Submitted Aug 2010  
}

\keyword{ package }

\seealso{
  \code{\link{fitCentipede}}, 
  \code{\link{plotProfile}}
}

\examples{
\dontrun{
#GETS EXAMPLE DATA FOR NRSF
data(NRSFcuts, package='CENTIPEDE')
data(NRSF_Anno, package='CENTIPEDE')

#FITS THE CENTIPEDE MODEL
centFit <- fitCentipede(Xlist = list(as.matrix(NRSFcuts)), Y=cbind(rep(1, dim(NRSF_Anno)[1]), NRSF_Anno[,5], NRSF_Anno[,6]))

#PLOTS IMAGE OF CUTSITES RANKED BY CENTIPEDE POSTERIORS
imageCutSites(NRSFcuts[order(centFit$PostPr),][c(1:100, (dim(NRSFcuts)[1]-100):(dim(NRSFcuts)[1])),])

#PLOT ESTIMATED FOOTPRINT
plotProfile(centFit$LambdaParList[[1]],Mlen=21)
}
}

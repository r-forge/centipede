\name{plotProfile} 
\alias{plotProfile} 
\title{Plots a footprint estimated by Centipede} 
\description{
  After fiting the CENTIPEDE model this function is useful for
  visualizing the footprint profile
}
\usage{
plotProfile(Profile, Mlen = 0, xlab = "Dist. to motif (bp)", ylab =
"Cut-site probability", legTitle = "", legPos = "topright")
}

\arguments{
  \item{Profile}{A vector with the profile estimated by CENTIPEDE}
  \item{Mlen}{Length of the motif for drawing vertical lines delimiting it}
  \item{xlab}{Label of the x axis}
  \item{ylab}{Label for the y axis}
  \item{legTitle}{Title for one of the plot corners}
  \item{legPos}{Corner where to place the legTitle}
}

\value{
 No return value, just the plot
}

\references{
  Roger Pique-Regi, Jacob F. Degner, Athma A. Pai, Daniel J. Gaffney,
  Yoav Gilad, Jonathan K. Pritchard. "Accurate inference of
  transcription factor binding from DNA sequence and chromatin
  accessibility data", Genome Research, Submitted Aug 2010  
}

\keyword{ package }

\seealso{
  \code{\link{fitCentipede}}
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


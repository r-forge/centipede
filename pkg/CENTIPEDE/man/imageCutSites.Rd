\name{imageCutSites} 
\alias{imageCutSites} 
\alias{imageCutSitesCombined} 
\title{Plots and image with the cleavage pattern} 
\description{
  After fiting the CENTIPEDE model this function is useful for
  visualizing the footprint cleavage pattern. 
}
\usage{
imageCutSites(data, main="", xlab = "Dist. to motif (bp)")
imageCutSitesCombined(data, main="",xlab= "Dist. to motif (bp)",breaks)
}

\arguments{
  \item{data}{A matrix with the data supplied to CENTIPEDE}
  \item{xlab}{Label of the x axis}
  \item{main}{Title for the plot}
  \item{breaks}{Breaks for the color assigments}
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

\keyword{ utils }

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


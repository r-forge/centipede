\name{fitCentipede} 
\alias{fitCentipede} 
\title{CENTIPEDE core function} 
\description{Fit the Centipede algorithm described in Pique-Regi, Degner, et al. (submitted).} 
\usage{ 
fitCentipede(Xlist,Y, LambdaParList, BetaLogit, NegBinParList, sweeps=50,DampLambda=0.0,DampNegBin=0.0,TrimP=0.0001,NRiter=5) 
} 
\arguments{ 
\item{Xlist}{A list of data-matrices representing discrete counts around candidate binding sites (e.g., DNaseI-seq, ChIP-seq for histone modifications, ...).  Each matrix in this list should have a row for each candidate binding site.  Numbers of columns may vary.} 
\item{Y}{A design matrix for logistic-prior component of CENTIPEDE model.  Each row corresponds to a candidate site.  Each column represents a score (e.g., conservation, distance to TSS) for that candidate site.  Use a column of 1s to fit an intercept term.} 
\item{LambdaParList}{An optional list matching the length of Xlist where starting values of Lambda parameters can be specified.} 
\item{BetaLogit}{An optional vector matching the columns of Y where starting values of Beta parameters of the logistic can be specified.} 
\item{NegBinParList}{An optional list matching the length of Xlist where starting values of each negative binomial distribution can be specified.}
\item{sweeps}{Maximum number of EM iterations to attempt.}
\item{DampLambda}{Level of shrinkage of Multinomial parameters (0-1).  0 represents no shrinkage and was generally found (for DNaseI data) to give highest specificity at the cost of some sensitivity.} 
\item{DampNegBin}{Level of shrinkage of Negative binomial parameters (0-1).}
\item{TrimP}{Artifically trims log likelihood ratios above this quantile of the total distribution.}
\item{NRiter}{Number of Newton-Raphson steps performed at each cycle of the EM for parameters without exact solutions.}
\item{\dots}{ not used. } 
}
\details{
We use a probabilistic framework 
known as a hierarchical mixture model that is described in greater detail in Pique-Regi et al (2010).  The likelihood for
a motif match l is written:

  P(D_l) = 
               P( Z_l=1| G_l ) * P(D_l|Z_l=1) + 
               P(Z_l=0|G_l) * P(D_l|Z_l=0)

where D_l and G_l represent
the observed experimental data and the prior information
around the motif match.
The data D_l are assumed to be generated from one of two
underlying distributions that form the mixture model. 
One distribution corresponds to the bound
state of transcription factors (Z_l=1) while the other distribution
corresponds to the unbound state (Z_l=0).

For each potential binding location `l', we calculate a prior
probability PI_l = P( Z_l=1 | G_l) that the site is bound by a TF. 
This prior probability is modeled using a logistic function:

  log(PI_l/(1-PI_l)) = Beta_0 + Beta_1 * (PWM Score)_l+ Beta_2 *(Cons. Score)_l + Beta_3 * (TSS Proximity)_l ... 

One could imagine many potential genomic annotations that might add additional information to this prior.

As experimental data ``D_l'', CENTIPEDE can combine multiple types of experiments 
The underlying assumption is that the different experiments can be considered 
conditionally independent given that the underlying state ``Z_l'' is known. 

For a given experimental data-type (e.g., DNase-seq),
the collection of reads in a region (200bp) around the motif
matches `l' can be represented by an L * S matrix
Each row corresponds to motif match location `l' and each column
`s' indexes the DNaseI cut position relative to the center and strand of this motif match.  
The total number of reads in the region is defined as R_l and is modeled with negative binomial distributions,
which depend on Alpha_1,Tau_1 for the bound class, and
Alpha_0,Tau_0 for the unbound class.  
While Poisson distributions may seem like the natural choice for the underlying process, 
the 2-parameter negative binomial distribution 
allows us to more accurately model the variance in sequence read rate. 
With these two
distributions we can capture open versus closed chromatin in DNaseI
hypersensitivity assays, or enrichment of certain histone
modifications associated with enhancers or repressors measured by
ChIP-seq assays. 
If the positional distribution of reads is not important 
(or not very informative) we can leave it unspecified  (i.e, any configuration is
equally likely).
This is the option we chose for the histone modification ChIP-seq
assays based on preliminary analysis showing that the read locations were
only weakly informative for these data (achieved by setting DampLambda=1). 
In contrast, for DNaseI the positional information can be very
informative as DNaseI leaves a distinctive
cleavage pattern (footprint) when Z_l=1. 
The spatial distribution of reads
surrounding the binding site is modeled with a multinomial
distribution
where the Lambda_l gives the probability that a read is obtained
from position index `s'.  
For Z_l=0, the TF is not bound,
so no specific footprint is expected. In this case
we find it works well to simply model the cut-site  distribution as uniform

The parameters of the CENTIPEDE model are estimated by maximizing the likelihood function 
using an expectation maximization (EM) algorithm (see paper for details). 
Once the model has converged, the posterior probability is used to infer
whether a TF is bound at location l. 

} 
\value{ 
A list containing elements:
\item{PostPr}{A vector of posterior probabilities that each site is in the bound class} 
\item{LogRatios}{A vector of the log likelihood  ratios between the bound and unbound models for each site - equivilent to the logistic transformation of the posterior probabilities (but with greater precision)} 
\item{PriorPr}{A vector of estimated prior probabilities that each site is in the bound class}
\item{PriorLogRatio}{A vector of the log prior ratios between the bound and unbound models for each site} 
\item{MultiNomLogRatio}{A vector of the log likelihood  ratios between the bound and unbound models coming from the multinomial component of the model only}
\item{NegBinLogRatio}{A vector of the log likelihood  ratios between the bound and unbound models coming from the negative binomial component of the model only}
\item{LambdaParList}{A vector of estimates for the Lambda parameters (e.g., the footprint profile)}
\item{BetaLogit}{A vector of estimated Beta coefficients of each element in the prior}
\item{LogLikEnd}{Final log-likelihood  of the data given the final model}
\item{NumIter}{The number of iterations EM algorithm took to converge}
} 
\author{Roger Pique-Regi and Jacob Degner} 
\examples{ 
#GETS EXAMPLE DATA FOR NRSF
data(NRSFcuts, package='Centipede2')
data(NRSF_Anno, package='Centipede2')

#FITS THE CENTIPEDE MODEL
centFit <- fitCentipede(Xlist = list(as.matrix(NRSFcuts)), Y=cbind(rep(1, dim(NRSF_Anno)[1]), NRSF_Anno[,5], NRSF_Anno[,6]))

#PLOTS IMAGE OF CUTSITES RANKED BY CENTIPEDE POSTERIORS
imageCutSites(NRSFcuts[order(centFit$PostPr),][c(1:100, (dim(NRSFcuts)[1]-100):(dim(NRSFcuts)[1])),])

#PLOT ESTIMATED FOOTPRINT
plotProfile(centFit$LambdaParList[[1]])

} 
\keyword{centipede} 


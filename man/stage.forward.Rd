\name{stage.forward}
\alias{stage.forward}
\title{
Nonlinear Forward stagewise regression using DCOL
}
\description{
The subroutine conducts forward stagewise regression using DCOL. Either DCOL roughening or spline roughening is conducted. 
}
\usage{
stage.forward(X, y, step.size = 0.01, stop.alpha = 0.01, 
stop.var.count = 20, roughening.method = "DCOL", tol = 1e-08, 
spline.df = 5, dcol.sel.only = FALSE, do.plot = F)
}
\arguments{
  \item{X}{
	The predictor matrix. Each row is a gene (predictor), each column is a sample. Notice the dimensionality is different than most other packages, where each column is a predictor. This is to conform to other functions in this package that handles gene expression type of data. 
}
  \item{y}{
	The numerical outcome vector.
}
  \item{step.size}{
	The step size of the roughening process. 
}
  \item{stop.alpha}{
	The alpha level (significance of the current selected predictor) to stop the iterations.	
}
  \item{stop.var.count}{
	The maximum number of predictors to select. Once this number is reached, the iteration stops.
}
  \item{roughening.method}{
	The method for roughening. The choices are "DCOL" or "spline". 
}
  \item{tol}{
	The tolerance level of sum of squared changes in the residuals. 	
}
  \item{spline.df}{
	The degree of freedom for the spline. 	
}
  \item{dcol.sel.only}{
	TRUE or FALSE. If FALSE, the selection of predictors will consider both linear and nonlinear association significance. 	
}
  \item{do.plot}{
	Whether to plot the points change in each step.	
}
}
\details{
	Please refer to the reference manuscript for details.
}
\value{
A list object is returned. The components include the following.
\item{found.pred}{ The selected predictors (row number).}
\item{ssx.rec}{The magnitude of variance explained using the current predictor at each step.}
\item{$sel.rec}{ The selected predictor at each step.}
\item{$p.rec}{ The p-value of the association between the current residual and the selected predictor at each step.}
}
\references{
	https://arxiv.org/abs/1601.05285
}
\author{
Tianwei Yu<tianwei.yu@emory.edu>
}

\seealso{
nvsd
}
\examples{
X<-matrix(rnorm(2000),ncol=20)
y<-sin(X[,1])+X[,2]^2+X[,3]
stage.forward(t(X),y,stop.alpha=0.001,step.size=0.05)
}
\keyword{ nonparametric }% use one of  RShowDoc("KEYWORDS")
\keyword{ variable selection }% __ONLY ONE__ keyword per line

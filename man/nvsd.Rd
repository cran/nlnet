\name{nvsd}
\alias{nvsd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Nonlinear Variable Selection based on DCOL	
}
\description{
This is a nonlinear variable selection procedure for generalized additive models. It's based on DCOL, using forward stagewise selection. In addition, a cross-validation is conducted to tune the stopping alpha level and finalize the variable selection.
}
\usage{
nvsd(X, y, fold = 10, step.size = 0.01, stop.alpha = 0.05, stop.var.count = 20, 
max.model.var.count = 10, roughening.method = "DCOL", do.plot = F, pred.method = "MARS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	The predictor matrix. Each row is a gene (predictor), each column is a sample. Notice the dimensionality is different than most other packages, where each column is a predictor. This is to conform to other functions in this package that handles gene expression type of data. 
}
  \item{y}{
	The numerical outcome vector.
}
  \item{fold}{
	The fold of cross-validation.
}
  \item{step.size}{
	The step size of the roughening process. 
}
  \item{stop.alpha}{
	The alpha level (significance of the current selected predictor) to stop the iterations.	
}
  \item{stop.var.count}{
	The maximum number of predictors to select in the forward stagewise selection. Once this number is reached, the iteration stops.
}
  \item{max.model.var.count}{
	The maximum number of predictors to select. Notice this can be smaller than the stop.var.count. Stop.var.count can be set more liniently, and this parameter controls the final maximum model size.
}
  \item{roughening.method}{
	The method for roughening. The choices are "DCOL" or "spline". 
}
  \item{do.plot}{
	Whether to plot the points change in each step.	
}
  \item{pred.method}{
	The prediction method for the cross validation variable selection. As forward stagewise procedure doesn't do prediction, a method has to be borrowed from existing packages. The choices include "MARS", "RF", and "SVM".
}
}
\details{
Please refer to the reference for details.
}
\value{
A list object is returned. The components include the following.
\item{selected.pred}{ The selected predictors (row number).}
\item{all.pred}{The selected predictors by the forward stagewise selection. The $selected.pred is a subset of this.}
}
\references{
	https://arxiv.org/abs/1601.05285
}
\author{
Tianwei Yu<tianwei.yu@emory.edu>
}

\seealso{
stage.forward
}
\examples{
X<-matrix(rnorm(2000),ncol=20)
y<-sin(X[,1])+X[,2]^2+X[,3]
nvsd(t(X),y,stop.alpha=0.001,step.size=0.05)
}
\keyword{ nonparametric }% use one of  RShowDoc("KEYWORDS")
\keyword{ variable selection }% __ONLY ONE__ keyword per line

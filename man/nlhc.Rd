\name{nlhc}
\alias{nlhc}
\title{
  Non-Linear Hierarchical Clustering
}
\description{
  The non-linear hierarchical clustering based on DCOL
}
\usage{
nlhc(array, hamil.method = "nn", concorde.path = NA, 
use.normal.approx = FALSE, normalization = "standardize", combine.linear = TRUE,
use.traditional.hclust = FALSE, method.traditional.hclust = "average")
}
\arguments{
  \item{array}{
the data matrix with no missing values
}
  \item{hamil.method}{
the method to find the hamiltonian path. 
}
  \item{concorde.path}{
If using the Concorde TSP Solver, the local directory of the solver
}
  \item{use.normal.approx}{
whether to use the normal approximation for the null hypothesis.
}
  \item{normalization}{
the normalization method for the array.
}
  \item{combine.linear}{
whether linear association should be found by correlation to combine with nonlinear association found by DCOL. 
}
  \item{use.traditional.hclust}{
whether traditional agglomerative clustering should be used. 
}
  \item{method.traditional.hclust}{
the method to pass on to hclust() if traditional method is chosen.
}
}
\details{
Hamil.method: It is passed onto the function tsp of library TSP. To use linkern method, the user needs to install concord as instructed in TSP.

use.normal.approx: If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.

normalization: There are three choices - "standardize" means removing the mean of each row and make the standard deviation one; "normal_score" means normal score transformation; "none" means do nothing. In that case we still assume some normalization has been done by the user such that each row has approximately mean 0 and sd 1.

combine.linear: The two pieces of information is combined at the start to initiate the distance matrix. 
}
\value{
Returns a hclust object same as the output of hclust().
Reference: help(hclust)
  \item{merge }{an n-1 by 2 matrix. Row i of merge describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm.}
  
  \item{height }{a set of n-1 real values, the value of the criterion associated with the clusterig method for the particular agglomeration}
  
  \item{order }{a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches.}
  
  \item{labels }{labels for each of the objects being clustered}
  \item{call }{the call which produced the result}
  \item{dist.method }{the distance that has been used to create d}
  \item{height.0}{original calculation of merging height}
}
\references{
http://www.ncbi.nlm.nih.gov/pubmed/24334400
}
\author{
 Tianwei Yu <tianwei.yu@emory.edu>
}

\seealso{
\code{\link{data.gen}}
}

\examples{
 ## generating the data matrix & hiden clusters as a sample
 input<-data.gen(n.genes=40, n.grps=4)
 ## now input includes data matrix and hiden clusters, so get the matrix as input.
 input<-input$data

 nlhc.data<-nlhc(input)
 plot(nlhc.data)
 ##get the merge from the input.
 merge<-nlhc.data$merge 
}
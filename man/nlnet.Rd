\name{nlnet}
\alias{nlnet}
\title{
  Non-Linear Network reconstruction from expression matrix
}
\description{
  Non-Linear Network reconstruction method
}
\usage{
nlnet(input, min.fdr.cutoff=0.05,max.fdr.cutoff=0.2, conn.proportion=0.007, 
gene.fdr.plot=FALSE, min.module.size=0, gene.community.method="multilevel", 
use.normal.approx=FALSE, normalization="standardize", plot.method="communitygraph")
}
\arguments{
  \item{input}{the data matrix with no missing values.
}
  \item{min.fdr.cutoff}{
  the minimun allowable value of the local false discovery cutoff in establishing links between genes. 
} \item{max.fdr.cutoff}{
   the maximun allowable value of the local false discovery cutoff in establishing links between genes.
} \item{conn.proportion}{
   the target proportion of connections between all pairs of genes, if allowed by the fdr cutoff limits. 
}
 \item{gene.fdr.plot}{
  whether plot a figure with estimated densities, distribution functions, and (local) false discovery rates.
}
  \item{min.module.size}{
  the min number of genes together as a module.
}
  \item{gene.community.method}{the method for community detection.
}
  \item{use.normal.approx}{whether to use the normal approximation for the null hypothesis. 
}
  \item{normalization}{the normalization method for the array. 
}
  \item{plot.method}{
  the method for graph and community ploting.
}
}
\details{
gene.community.method: It provides three kinds of community detection method:
  "mutilevel", "label.propagation" and "leading.eigenvector".
  
use.normal.approx: If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.

normalization: There are three choices:
"standardize" means removing the mean of each row and make the standard deviation one; 
"normal_score" means normal score transformation; 
"none" means do nothing. In that case we still assume some normalization has been done by the user such that each row has approximately mean 0 and sd 1. 

plot.method: It provides three kinds of ploting method:
  "none" means ploting no graph,
  "communitygraph" means ploting community with graph, 
  "graph" means ploting graph, 
  "membership" means ploting membership of the community 
}
\value{
  it returns a graph and the community membership of the graph.
 \item{algorithm}{The algorithm name for community detection}
 \item{graph}{An igraph object including edges : Numeric vector defining the edges, the first edge points from the first element to the second, the second edge from the third to the fourth, etc.}
 \item{community}{Numeric vector, one value for each vertex, the membership vector of the community structure.}
}
\author{
  Haodong Liu <liuhaodong0828@gmail.com>
}


\seealso{
  \code{\link{data.gen}}
}
\examples{
 
 ## generating the data matrix & hiden clusters as a sample
 input<-data.gen()
 ##Now input includes data matrix and hiden clusters, so get the matrix as input.
 input<-input$data
 
 ##change the ploting method
 nlnet(input,plot.method="graph")
 ## see fdrtool image
 nlnet(input,gene.fdr.plot=TRUE)
 
 ## use different community detection method
 nlnet(input,gene.community.method="label.propagation")
 
 ## change the fdr pro to control connections of genes
 ## adjust the modularity size
 nlnet(input,conn.proportion=0.005,min.module.size=10)
 
 ## get the result and see it values
 result<-nlnet(input)
 graph<-result$graph ##a igraph object.
 comm<-result$community ##community of the graph
}

\name{KPC}
\alias{KPC}
\title{
 implementation of K-Profiles Clustering
}
\description{
 implementation of K-Profiles Clustering
}
\usage{
 KPC(dataset, nCluster, maxIter = 100, p.max = 0.2, p.min = 0.05)
}
\arguments{
  \item{dataset}{
the data matrix with genes in the row and samples in the column
}
  \item{nCluster}{
the number of clusters K
}
  \item{maxIter}{
the maximum number of iterations
}
  \item{p.max}{
the starting p-value cutoff to exclude noise genes
}
  \item{p.min}{
the final p-value cutoff to exclude noise genes
}
}
\value{
 Return a list about gene cluster and the list of value p
 \item{cluster}{gene cluster}
 \item{p.list}{a list of value p}
}
\references{
http://www.hindawi.com/journals/bmri/aa/918954/
}
\author{
Tianwei Yu <tianwei.yu@emory.edu>
}
\seealso{
 \code{\link{data.gen}}
}
\examples{

 ## generating the data matrix & hiden clusters as a sample
 input<-data.gen()
 ## now input includes data matrix and hiden clusters, so get the matrix as input.
 input<-input$data
 
 ## set nCluster value to 5
 kpc<-KPC(input,nCluster=5)
 
 ## set maximun number of iteration to 200
 kpc<-KPC(input,nCluster=10,maxIter=200)
 
 ##get the hiden cluster result from "KPC"
 cluster<-kpc$cluster
 ##get the list of p
 p<-kpc$p.list
}
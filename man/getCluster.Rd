\name{getCluster}
\alias{getCluster}
\title{Detect clusters of genomic features (i.e TFBS) in a defined window size
following user defined options}
\usage{
getCluster(x, w, c, overlap = 0, greedy = TRUE, chr = NULL, s = ".",
order = NULL, n_cores = 2)

}
\description{
Given a data frame or bed/vcf files, this function will look for
clusters of genomic features found within a user defined window size and
satisfiying user defined categorical conditions.
}

\arguments{
    \item{x}{Vector of file names or a data frame}

    \item{w}{Value for the desired cluster size (numeric)}

    \item{c}{A named vector for condition relative to each site name}

    \item{overlap}{numeric; If negative, it will correspond to the maximum
                overlap allowed between consecutive clusters, if positive,
                it will correspond to the minimum gap allowed between
                consecutive clusters. The default is set to 0}


    \item{greedy}{logical; If FALSE, the formed clusters will
            contain the exact combination of sites defined in "c".
            If TRUE, clusters containing the condition and more will
            be labeled as TRUE.The default is set to TRUE}

    \item{chr}{character; chromosome to cluster on. The default is set to
                NULL so the clustering will be performed on all chromosomes}

    \item{s}{character; Strand to cluster on. Can be "+", "-" or
                ".". The default is set to "."}

    \item{order}{vector; A vector containing the order of sites the user
                        wants the clusters to contains. By default it is
                        set to NULL}

    \item{n_cores}{numeric; The number of cores to be used in parallel.
                    By default it is set to 2 cores}



}

\details{
The function \emph{getCluster} will cluster coordinates based on a user defined
number of sites and window size.
The user needs to specify the condition for clustering, which is the
sites and the number of sites required in each cluster. The user can
also set the distance required between consective clusters using the
overlap argument. If overlap is a negative, this will represent the
maximum overlap allowed between clusters. If overlap is positive, this
means that the clusters should have a minimum gap of the given value.

The user can also choose to cluster on a specific strand by using the strand
option, and specific chromosome by passing on the chromosome name as an
argument.

\emph{x} is a vector of files or a data frame.
In the case of the data frame it should have 5 columns
\emph{chr start end strand  site}. \emph{start} and \emph{end} column
need to be numeric or integer, while the rest of the columns are of type
character.

The Greedy option controls if more sites are allowed in the
cluster. When set to TRUE, getCluster will return sites containing
the required number of sites and more as long as the window size
condition is satisfied. When set to FALSE, getCluster will return
clusters containing the exact number of the required sites. Default is
set to TRUE.

\emph{order} allows the user to choose whether the clusters should contain
sites ordered in a certain way. When set to NULL, sites order is not
important.
When defined, the order of the sites in the cluster must abide by the
user defined order. Default is NULL.
}

\examples{
 x = data.frame(chr = c("chr1","chr2","chr1","chr1","chr2","chr2"),
 start = c(2,5,7,10,12,14),end = c(5,10,12,14,18,20),
 strand = c("+","+","+","+","+","+" ),
 site = c("s1", "s2", "s1","s2","s1", "s2"))

clusters = getCluster(x, w = 10, c = c("s1"=1,"s2"=2),
greedy = TRUE, overlap = -5, s = "+", order = c("s1","s2","s1"))

}


\value{
The returned value is a data frame that contains the clusters found.
}

\author{
  Pierre Khoueiry

  Ghiwa Khalil
}

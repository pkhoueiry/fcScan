\name{getCluster}
\alias{getCluster}
\title{Detect clusters of genomic features (i.e TFBS) in a specified window size
and using user defined options}
\usage{
getCluster(x, w, c, overlap = 0, greedy = TRUE, seqnames = NULL, s = ".",
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
            be labeled as TRUE as long as window size is below or equal
            the defined value. The default is set to TRUE}

    \item{seqnames}{character; seqname/chromosome to cluster on. The default is set to
                NULL so the clustering will be performed on all seqnames}

    \item{s}{character; Strand to cluster on. Can be "+", "-" or
                ".". The default is set to "." which indicate both clusters}

    \item{order}{vector; A vector containing the order of sites the user
                        wants the clusters to contain. By default it is
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
option, and specific seqname by specifying the seqnames as an
argument.

\emph{x} is a vector of files or a data frame.
When data frame is given as input, it should have 5 columns
\emph{seqnames start end strand site}. \emph{start} and \emph{end} column
need to be numeric or integer and the rest of the columns are of type
character.

The Greedy option controls if more sites are allowed in the
cluster. When set to TRUE, getCluster will return clusters containing
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
x = data.frame(seqnames = rep("chr1", times = 16),
start = c(10,17,25,27,32,41,47,60,70,87,94,99,107,113,121,132),
end = c(15,20,30,35,40,48,55,68,75,93,100,105,113,120,130,135),
strand = rep("+", 16),
site = c("s1","s2","s2","s1","s2","s1","s1","s2","s1","s2","s2","s1","s2","s1","s1","s2"))

clusters = getCluster(x, w = 25, c = c("s1"=1,"s2"=2),
greedy = TRUE, overlap = -5, s = "+", order = c("s1","s2","s1"))

}


\value{
The returned value is a GRanges object containing the clusters.
}

\author{
  Pierre Khoueiry

  Ghiwa Khalil
}

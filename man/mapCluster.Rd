\name{mapCluster}
\alias{mapCluster}
\title{Identify orthologous coordinates between species}
\usage{
mapCluster(x, kent, chrom.sizes, chain)
}
\description{
Detect orthologous regions between species using pslMap from kent utilities.
}

\arguments{
    \item{x}{data frame of BED format containing the clusters to be mapped}

    \item{kent}{path to kentUtils. NULL by default}

    \item{chrom.sizes}{Path to chrom.sizes file}

    \item{chain}{Path to chain file}

}

\details{
The function \emph{mapCluster} uses bedToPsl and pslMap from kentUtils to detect
orthologous regions between 2 species.kentUtils are assumed to in the
system path. Otherwise, the path needs to be specified.
}

\examples{
##clusters is a dataframe that contains clusters to be mapped
clusters = data.frame(chr = c("chr22", "chr22", "chr22", "chr22", "chr22"
), start = c(22820958, 31964079, 36411953, 44469038, 49453465
), end = c(22826500, 31981177, 36426581, 44483647, 49463799),
    id = c("c89", "c214", "c270", "c367", "c438"), score = c(1,
    1, 1, 1, 1), strand = c("+", "+", "+", "+", "+"))

mapCluster( x = clusters, kent = NULL,
chrom.sizes = "gfScan/extdata/hg19.chrom.sizes",
chain = "gfScan/extdata/hg19ToMm10Chr22.over.chain")
}

\value{
The returned value is a data frame that contains the corresponding mapped region
for each cluster.
}

\author{
Pierre Khoueiry
Ghiwa Khalil
}

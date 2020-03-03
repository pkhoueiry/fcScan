fcScan
======

fcScan is a Bioconductor package aiming at clustering genomic features based on a defined window size and combination of sites.

Dependencies
-----------

The below R packages are required for installation:

+ stats
+ plyr
+ utils
+ VariantAnnotation
+ SummarizedExperiment
+ rtracklayer
+ GenomicRanges
+ IRanges

Installation
------------

1. Git clone the project directory `git clone https://github.com/pkhoueiry/fcScan.git`
2. From the terminal run `R CMD build fcScan`
3. Run `R CMD check fcScan_1.0.1.tar.gz`
4. Run  `R CMD BiocCheck fcScan_1.0.1.tar.gz`
5. Install package by `R CMD INSTALL fcScan_1.0.1.tar.gz`
6. Start R and load library using `library(fcScan)`




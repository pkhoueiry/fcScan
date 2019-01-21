fcScan
======

fcScan is a Bioconductor package aiming at clustering genomic features based on a defined window size and combination of sites.

Dependencies
-----------

The below R packages are required for installation:

+ parallel
+ doParallel
+ foreach
+ stats
+ utils
+ VariantAnnotation
+ rtracklayer
+ matrixStats
+ GenomicRanges
+ tools


Installation
------------

1. Git clone the project directory `git clone https://gitlab.com/pklab/fcScan.git`
2. From the terminal run `R CMD build fcScan`
3. Run `R CMD check fcScan_0.99.0.tar.gz`
4. Run  `R CMD BiocCheck fcScan_0.99.0.tar.gz`
5. Intsall package by `R CMD INSTALL fcScan_0.99.0.tar.gz`
6. Switch to R and load library by `library(fcScan)`




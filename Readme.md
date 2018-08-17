## All the files required for building the package are found in the master branch of SECOMOD.


1. Clone git repository these files into a directory
2. From the terminal run `R CMD build fcScan`
3. Run `R CMD check fcScan_0.99.0.tar.gz`
4. Run  `R CMD BiocCheck fcScan_0.99.0.tar.gz`
5. Switch to R and load library by `library(fcScan)`
6. Intsall package by `R CMD INSTALL fcScan_0.99.0.tar.gz`

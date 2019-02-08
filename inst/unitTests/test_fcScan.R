test_fcScan <- function() {
    clusters <- GRanges(seqnames = c("chr1","chr1","chr1","chr1","chr1","chr1","chr1"
                                     ,"chr1","chr1","chr1","chr1","chr1","chr1"), 
                  ranges= IRanges(start= c(11,18,26,48,61,71,88,95,100,108,114,122,133),
                                  width = c(25,23,23,21,15,23,18,19,21,23,22,14,3)),
                  size= c(25,23,23,21,15,23,18,19,21,23,22,14,3),
                  isCluster = c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,
                                FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                  id = c("c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13"),
                  status = c("orderFail","orderFail","PASS","combnFail","combnFail","combnFail",
                             "orderFail","orderFail","combnFail","combnFail","combnFail","combnFail",
                             "combnFail"),
                sites = c("s1,s2,s2,s1" ,"s2,s2,s1,s2" ,"s2,s1,s2,s1", "s1,s2","s2,s1",
                          "s1,s2","s2,s2,s1","s2,s1,s2","s1,s2,s1","s2,s1,s1",
                          "s1,s1,s2","s1,s2","s2"),
                score =c(1,1,1,1,1,1,1,1,1,1,1,1,1)
)
    
checkEquals(getCluster(x, w = 25, c = c("s1"=1,"s2"=2),
            greedy = TRUE, overlap = -5, s = "+", order = c("s1","s2","s1"))
            , clusters)
}


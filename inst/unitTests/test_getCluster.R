require(GenomicRanges)
x1 = data.frame(seqnames = rep("chr1", times = 16),
               start = c(10,17,25,27,32,41,47,60,70,87,94,99,107,113,121,132),
               end = c(15,20,30,35,40,48,55,68,75,93,100,105,113,120,130,135),
               strand = rep("+", 16),
               site = c("s1","s2","s2","s1","s2","s1","s1","s2","s1","s2","s2","s1","s2","s1","s1","s2"))

t1 = GRanges(
    seqnames = Rle("chr1", 13),
    ranges = IRanges(c(11L, 18L, 26L, 48L, 61L, 71L, 88L, 95L, 100L, 108L, 114L, 122L, 133L),
                     end = c(35L, 40L, 48L, 68L, 75L, 93L, 105L, 113L, 120L, 130L, 135L, 135L, 135L)),
    strand = Rle("+", 13),
    sites = as.character(c("s1,s2,s2,s1", "s2,s2,s1,s2", "s2,s1,s2,s1", "s1,s2", "s2,s1", "s1,s2", "s2,s2,s1",
             "s2,s1,s2", "s1,s2,s1", "s2,s1,s1", "s1,s1,s2", "s1,s2", "s2")),
    isCluster = as.logical(Rle(c(FALSE, TRUE, FALSE), c(2, 1, 10))),
    status = as.character(Rle(c("orderFail", "PASS", "combnFail", "orderFail", "combnFail"), c(2, 1, 3, 2, 5)))
)

## testing cases with one cluster as output
x2 = data.frame(seqnames = rep("chr1", times = 1),
               start = c(10,17,25),
               end = c(15,20,30),
               strand = rep("+", 3),
               site = c("s1","s2","s2"))

t2 = GRanges(
    seqnames = Rle("chr1", 1),
    ranges = IRanges(c(11L), end = c(30L)),
    strand = Rle("+", 1),
    sites = as.character(c("s1,s2,s2")),
    isCluster = as.logical(Rle(TRUE, 1)),
    status = as.character(Rle(c("PASS"), 1)))


## tests
test_getCluster <- function() {
    ## Generic case
    checkEquals(getCluster(x1, w = 25, c = c("s1"=1,"s2"=2), greedy = TRUE, overlap = -5, s = "+", order = c("s1","s2","s1"), verbose = TRUE), t1)
    ## One cluster found
    checkEquals(getCluster(x2, w = 20, c = c("s1"=1,"s2"=2), greedy = TRUE, overlap = 0, s = "+", order = c("s1","s2","s2"), verbose = TRUE), t2)
    ## Zero cluster found
    checkEquals(getCluster(x2, w = 15, c = c("s1"=1,"s2"=2), greedy = TRUE, overlap = 0, s = "+", order = c("s1","s2","s2"), verbose = FALSE), NULL)
}




#library(GenomicRanges)
#library(RUnit)
require(GenomicRanges)

x1 = data.frame(seqnames = rep("chr1", times = 17),
               start = c(1,10,17,25,27,32,41,47,60,70,87,94,99,107,113,121,132),
               end = c(8,15,20,30,35,40,48,55,68,75,93,100,105,113,120,130,135),
               strand = rep("+", 17),
               site = c("s3","s1","s2","s2","s1","s2","s1","s1","s2","s1","s2","s2","s1","s2","s1","s1","s2"))

t1 = GRanges(
    seqnames = Rle("chr1", 3),
    ranges = IRanges(c(11L, 26L, 95L),
        end = c(20L, 35L, 105L)),
    strand = Rle("+", 3),
    sites = as.character(c("s1,s2", "s2,s1", "s2,s1")),
    isCluster = as.logical(c(TRUE, FALSE, FALSE)),
    status = as.character(c("PASS", "orderFail", "orderFail"))
)

x2 = data.frame(seqnames = rep("chr1", times = 3),
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
    status = as.character(Rle(c("PASS"), 1))
)

t3 = GRanges(
    seqnames = Rle("chr1", 3),
    ranges = IRanges(c(11L, 26L, 95L), end = c(20L, 35L, 105L)),
    strand = Rle("+", 3),
    sites = as.character(c("s1,s2", "s2,s1", "s2,s1")),
    isCluster = as.logical(Rle(TRUE, 3)),
    status = as.character(Rle(c("PASS"), 3))
)

t4 = GRanges(
    seqnames = Rle("chr1", 1),
    ranges = IRanges(c(2L), end = c(15L)),
    strand = Rle("+", 1),
    sites = as.character(c("s3,s1")),
    isCluster = as.logical(Rle(TRUE, 1)),
    status = as.character(Rle(c("PASS"), 1))
)


t6 = GRanges(
    seqnames = Rle("chr1",2),
    ranges = IRanges(c(11L, 61L), end = c(55L, 105L)),
    strand = Rle("+", 2),
    sites = as.character(c("s1,s2,s2,s1,s2,s1,s1","s2,s1,s2,s2,s1")),
    isCluster = as.logical(c(TRUE, TRUE)),
    status = as.character(c("PASS", "PASS"))
)

t7 = GRanges(
    seqnames = Rle("chr1",3),
    ranges = IRanges(c(11L, 61L, 71L), end = c(40L, 105L, 113L)),
    strand = Rle("+", 3),
    sites = as.character(c("s1,s2,s2,s1,s2","s2,s1,s2,s2,s1","s1,s2,s2,s1,s2")),
    isCluster = as.logical(c(TRUE, FALSE, TRUE)),
    status = as.character(c("PASS", "orderFail", "PASS"))
)

t8 = GRanges(
    seqnames = Rle("chr1", 4),
    ranges = IRanges(c(11L, 18L, 71L, 88L), end = c(55L, 68L, 120L, 135L)),
    strand = Rle("+",4),
    sites = as.character(c("s1,s2,s2,s1,s2,s1,s1", "s2,s2,s1,s2,s1,s1,s2", "s1,s2,s2,s1,s2,s1", "s2,s2,s1,s2,s1,s1,s2")),
    isCluster = as.logical(c(FALSE, TRUE, FALSE, TRUE)),
    status = as.character(c("orderFail", "PASS", "orderFail", "PASS"))
)

t9 = GRanges(
    seqnames = Rle("chr1", 6),
    ranges = IRanges(c(11L, 26L, 61L, 95L, 108L, 122L), end = c(20L, 35L, 75L, 105L, 120L, 135L)),
    strand = Rle("*", 6),
    sites = as.character(c("s1,s2", "s2,s1", "s2,s1", "s2,s1", "s2,s1", "s1,s2")),
    isCluster = as.logical(Rle(TRUE, 6)),
    status = as.character(Rle("PASS", 6))
)

t10 = GRanges(
    seqnames = Rle("chr1", 7),
    ranges = IRanges(c(11L, 26L, 33L, 61L, 95L, 108L, 122L), end = c(20L, 35L, 48L, 75L, 105L, 120L, 135L)),
    strand = Rle("*", 7),
    sites = as.character(c("s1,s2", "s2,s1", "s2,s1", "s2,s1", "s2,s1", "s2,s1", "s1,s2")),
    isCluster = as.logical(Rle(TRUE, 7)),
    status = as.character(Rle("PASS", 7))
)

t11 = GRanges(
    seqnames = Rle("chr1", 3),
    ranges = IRanges(c(11L, 26L, 95L), end = c(20L, 35L, 105L)),
    strand = Rle("*", 3),
    sites = as.character(c("s1,s2", "s2,s1", "s2,s1")),
    isCluster = as.logical(Rle(TRUE, 3)),
    status = as.character(Rle("PASS", 3))
)

x3 = data.frame(seqnames = rep("chr1", times = 6),
               start = c(10,17,25,40,42,55),
               end = c(15,20,30,50,56,70),
               strand = c("+", "-", "-", "-", "+", "+"),
               site = c("s1","s2","s2","s2","s1","s1"))

t12 = GRanges(
    seqnames = Rle("chr1", 2),
    ranges = IRanges(c(11L, 41L), end = c(20L, 56L)),
    strand = Rle("*", 2),
    sites = as.character(c("s1,s2", "s2,s1")),
    isCluster = as.logical(Rle(TRUE, 2)),
    status = as.character(Rle("PASS", 2))
)

x4 = data.frame(seqnames = rep("chr1", times = 6),
               start = c(10,17,25,40,42,55),
               end = c(15,20,30,50,56,70),
               strand = c("+", "+", "+", "+", "+", "+"),
               site = c("s1","s2","s3","s2","s1","s3"))

t13 = GRanges(
    seqnames = Rle("chr1", 2),
    ranges = IRanges(c(11L, 41L), end = c(20L, 56L)),
    strand = Rle("*", 2),
    sites = as.character(c("s1,s2", "s2,s1")),
    isCluster = as.logical(Rle(TRUE, 2)),
    status = as.character(Rle("PASS", 2))
)

t14 = GRanges(
    seqnames = Rle("chr1", 1),
    ranges = IRanges(c(11L), end = c(30L)),
    strand = Rle("*", 1),
    sites = as.character(c("s1,s2,s3")),
    isCluster = as.logical(Rle(FALSE, 1)),
    status = as.character(Rle("ExcludedSites", 1))
)

t15 = GRanges(
    seqnames = Rle("chr1", 2),
    ranges = IRanges(c(11L, 18L), end = c(30L, 30L)),
    strand = Rle("*", 2),
    sites = as.character(c("s1,s2,s3", "s2,s3")),
    isCluster = as.logical(c(FALSE,TRUE)),
    status = as.character(c("ExcludedSites", "PASS"))
)

x5 = data.frame(seqnames = rep("chr1", times = 10),
               start = c(10,17,25,40,42,55,72,75,90,95),
               end = c(15,20,30,50,56,70,80,89,100,115),
               strand = c("+", "+", "+", "+", "+", "+", "+", "+", "+", "+"),
               site = c("s1","s2","s3","s2","s1","s3","s1","s3","s2","s1"))

t16 = GRanges(
    seqnames = Rle("chr1", 2),
    ranges = IRanges(c(11L, 18L), end = c(50L, 56L)),
    strand = Rle("*", 2),
    sites = as.character(c("s1,s2,s3,s2", "s2,s3,s2,s1")),
    isCluster = as.logical(c(FALSE, FALSE)),
    status = as.character(c("ExcludedSites", "ExcludedSites"))
)

x6 = data.frame(seqnames = rep("chr1", times = 10),
               start = c(10,17,25,40,42,55,72,75,90,95),
               end = c(15,20,30,50,56,70,80,89,100,112),
               strand = c("+", "+", "+", "+", "+", "+", "+", "+", "+", "+"),
               site = c("s1","s2","s2","s3","s1","s3","s1","s2","s2","s1"))

t17 = GRanges(
    seqnames = Rle("chr1", 3),
    ranges = IRanges(c(11L, 18L, 73L), end = c(50L, 56L, 112L)),
    strand = Rle("*", 3),
    sites = as.character(c("s1,s2,s2,s3", "s2,s2,s3,s1", "s1,s2,s2,s1")),
    isCluster = as.logical(c(FALSE, FALSE, TRUE)),
    status = as.character(c("ExcludedSites", "ExcludedSites", "PASS"))
)

x7 = GRanges(
    seqnames = Rle("chr1", 16),
    ranges = IRanges(start = c(10L,17L,25L,27L,32L,41L,47L,60L,70L,87L,94L,99L,107L,113L,121L,132L),
    end = c(15L,20L,30L,35L,40L,48L,55L,68L,75L,93L,100L,105L,113L,120L,130L,135L)),
    strand = Rle("+",16),
    site = c("s1","s2","s2","s1","s2","s1","s1","s2",
        "s1","s2","s2","s1","s2","s1","s1","s2"))

t18 = GRanges(
    seqnames = Rle("chr1", 2),
    ranges = IRanges(start = c(10L, 87L), end = c(30L, 105L)),
    strand = Rle("*",2),
    sites = as.character(c("s1,s2,s2", "s2,s2,s1")),
    isCluster = as.logical(c(TRUE, TRUE)),
    status = as.character("PASS","PASS")
)

test_getCluster <- function() {
    #general test t1#
    checkEquals(getCluster(x1, w = 11, c = c("s1" = 1, "s2" = 1), greedy = FALSE, order = c("s1","s2"), s = "+", verbose = TRUE), t1)
    #same previous test t1, but specifying sites "-" - so expecting NULL#
    checkEquals(getCluster(x1, w = 11, c = c("s1" = 1, "s2" = 1), greedy = FALSE, order = c("s1","s2"), s = "-", verbose = TRUE), NULL)
    #test t2 - get ONE cluster only#
    checkEquals(getCluster(x2, w = 20, c = c("s1" = 1,"s2" = 2), greedy = TRUE, overlap = 0, s = "+", order = c("s1","s2","s2"), verbose = TRUE), t2)
    #test t3 - get TWO clusters only#
    checkEquals(getCluster(x1, w = 10, c = c("s1" = 1, "s2" = 1), greedy = FALSE, s = "+", verbose = TRUE), t3)
    #test t4 - get one cluster having only sites s1 and s3#
    checkEquals(getCluster(x1, w = 15, c = c("s1" = 1, "s3" = 1), greedy = TRUE, s = "+", verbose = TRUE), t4)
    #test t5 - get NULL#
    checkEquals(getCluster(x1, w = 10, c = c("s1" = 1, "s3" = 1), greedy = TRUE, order = c("s1", "s3"), s = "+", verbose = TRUE), NULL)
    #test t6 - get only two clusters and zero orderFail#
    checkEquals(getCluster(x1, w = 50, c = c("s1" = 2, "s2" = 3), greedy = TRUE, s = "+", verbose = FALSE), t6)
    #test t7 - get two clusters and one orderFail#
    checkEquals(getCluster(x1, w = 50, c = c("s1" = 2, "s2" = 3), greedy = FALSE, order = c("s1","s2","s2","s1","s2"), s = "+", verbose = TRUE), t7)
    #test t8 - get two clusters and two orderFail#
    checkEquals(getCluster(x1, w = 50, c = c("s1" = 2, "s2" = 2), greedy = TRUE, order = c("s2","s1","s1","s2"), s = "+", verbose = TRUE), t8)
    #test t9 - get clusters not overlapping#
    checkEquals(getCluster(x1, w = 16, c = c("s1" = 1, "s2" = 1), greedy = FALSE, verbose = TRUE), t9)
    #test t10 - get overlapping clusters#
    checkEquals(getCluster(x1, w = 16, c = c("s1" = 1, "s2" = 1), greedy = FALSE, overlap = -3, verbose = TRUE), t10)
    #test t11 - get clusters with minimum gap of 6#
    checkEquals(getCluster(x1, w = 10, c = c("s1" = 1, "s2" = 1), greedy = FALSE, overlap = 6, verbose = TRUE), t11)
    #test t12 - get two clusters regardless of strand information#
    checkEquals(getCluster(x3, w = 16, c = c("s1" = 1, "s2" = 1), greedy = TRUE, overlap = -2, verbose = TRUE), t12)
    #test t13 - regardless of strand information#
    checkEquals(getCluster(x3, w = 16, c = c("s1" = 1, "s2" = 1), greedy = TRUE, verbose = TRUE), t13)
    #test t14 - get clusters of s1 and s3 sites with s2 excluded#
    checkEquals(getCluster(x4, w = 20, c = c("s1" = 1, "s2" = 0, "s3" = 1), greedy = TRUE, verbose = TRUE), t14)
    #test t15 - get one TRUE cluster and one with excluded sites#
    checkEquals(getCluster(x4, w = 20, c = c("s1" = 0, "s2" = 1, "s3" = 1), greedy = TRUE, verbose = TRUE), t15)
    #test t16 - get two ExcludedSites#
    checkEquals(getCluster(x5, w = 40, c = c("s1" = 1, "s2" = 2, "s3" = 0), greedy = TRUE, verbose = TRUE), t16)
    #test t17 - similar to t16 but Greedy = FALSE
    checkEquals(getCluster(x5, w = 40, c = c("s1" = 1, "s2" = 2, "s3" = 0), greedy = FALSE, verbose = TRUE), NULL)
    #test t18 - get three consecutive ExcludedSites and One TRUE cluster#
    checkEquals(getCluster(x6, w = 40, c = c("s1" = 1, "s2" = 2, "s3" = 0), order = c("s1", "s2", "s2"), greedy = TRUE, verbose = TRUE), t17)
    #test t19 - input GRanges#
    checkEquals(getCluster(x7, w = 25, c = c("s1"=1,"s2"=2)), t18)
}

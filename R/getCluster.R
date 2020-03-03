##' @title fcScan

##' @export getCluster




getCluster <- function(x, w, c, overlap = 0, greedy = FALSE, seqnames = NULL,
                    s = "*" , order = NULL, sites_orientation = NULL, 
                    intra_distance = 0, verbose = FALSE) {

    sitesToExclude <- NULL
    final = NULL
    start.time = Sys.time()

    ##check arguments
    if (missing(x)) {
        stop("A data frame, Granges object or input files are needed")
    }

    ## Checking input format
    if (is(x, "data.frame") || is(x, "GRanges")) {
        ## print("data.frame input")
        if(is.null(names(c))) {
            stop("When input is a data frame or GRanges object, site names
in condition must be explicitly defined")
        }
    }
        
        if(is.data.frame(x)){
        x = makeGRangesFromDataFrame(x, keep.extra.columns = TRUE,
            starts.in.df.are.0based = TRUE)
    }

    else if(is(x, "GRanges")){
        x
    }    
      
    else {
        ## print("File(s) input")
        if (length(c) != length(x))
            stop("Condition and files should be of same length")
        
        if(length(which(names(c)=="")>0)){
            stop("You either give names to all conditions or leave them empty")
        }

        ##assigning names to condition if null
        if (is.null(names(c))) {
            names(c) = seq(from = 1,to = length(c),by = 1)
        }

        x = load_data(all_files = x, c = c)
        x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE,
            starts.in.df.are.0based = TRUE)
    }
    
    if(length(which(duplicated(names(c))))!=0){
        stop("Names of conditions must be unique")
    }

    if(length(which(c<0))!=0 | !(is.numeric(c))){
        stop("Only positive integers allowed")
    }

    ## overlap will accept integers only
    if(!(overlap%%1 == 0)){
        stop("Only integers are allowed for overlap")
    }

    indexToExclude <- which(c == 0)

    ## Getting conditions to exclude (i.e with 0 values)
    if(length(indexToExclude) != 0){
        sitesToExclude <- names(indexToExclude)
    }
    
    ## Test if window is a positive integer
    if (w < 0 | w%%1!=0)
        stop("Window size should be a positive integer")
    
    if (!(s %in% c("+", "-", "*")))
        stop('Strand needs to be a valid option. Accepted options are "+", "-"
                and "*"')

    ##if not greedy and order is given
    if(is.null(order) == FALSE){
    ##checking for consistency of name in order and condition
    if(!(all(order %in% names(c))))
        stop("site names between order and condition do not match")

    ##check number of sites if equal in condition and order
    if(greedy == FALSE){
        count_elements <- c(count(order)$freq)
        names(count_elements) <- count(order)$x
        if(!(all(sort(count_elements) == sort(c)))){
            stop("Greedy is FALSE and order is larger than condition")
            }
        }
    }
    
    #site orientation works only if order is provided
    if(!(is.null(sites_orientation)) & is.null(order)){
        stop("sites_orientation cannot be used unless order is specified")
    }

    #site orientation should be positive or negative
    if(!(all(sites_orientation %in% c("+","-")))){
        stop('Strand should be only "+" or "-"')
    }

    #site orientation input length should have same order length
    if(length(order) != length(sites_orientation) & 
        !(is.null(sites_orientation))){
        stop("Orientation must be added to all sites in 'order' respectively")
    }

    ##intra_distance should be either positive, negative or zero
    if(!(is.numeric(intra_distance) || intra_distance == 0)){
        stop("Only integers allowed")
    }

    ##check verbose input argument
    if( !(verbose %in% c("TRUE", "FALSE"))) {
        stop("verbose should be TRUE or FALSE")
    }

    n = sum(c)

    ##getting sites found on the required seqnamesom
    ### change 2: subsetting from GRanges according to seqnames ###
    if (length(seqnames) > 0) {
        x = x[seqnames(x) %in% seqnames]
    }

    ##getting sites found on the required strand
    ### change 3: subsetting from GRanges according to strand ###
    if (s != "*") {
        x <- x[strand(x) == s]
    }

    ##check if the sites given in condition c are found in the data'
    if( !all(names(c) %in% x$site)) {
        message("Sites in condition do not match sites in data")
        return(NULL)
    }

    cat (length(x), " entries loaded", "\n")

    ##need to subset to keep only the sites required by the user,
    ##if the user wants sites "a","b" and input has "a", "b" and "c",
    ##we subset to get sites "a", "b"
    
    if(!(is.null(sitesToExclude))){
        x
    }
    
    else{
    x = subset(x, x$site %in% names(c))
    }
    
    unique_seqnames = unique((seqnames(x)))

    ## creating an array to fill results
    res = array(data = NA, dim = c(length(x), ncol = 7))
    colnames(res) = c("seqnames", "start", "end", "site", "strand",
                        "isCluster", "status")


    ## looping over chromosomes 
    for(seq in seq_along(unique_seqnames)){

        gr = subset(x, seqnames == unique_seqnames[seq])
        gr = sort(gr, ignore.strand = TRUE)

        if (length(gr) >= n) {
            result = cluster_sites(gr, w, c, overlap, n,
                    res, s, greedy, order, sitesToExclude, 
                        sites_orientation, intra_distance)
            final = rbind(final, result)
        }
    }

    if(length(final) !=0 ){
        final <- GRanges(
            seqnames = as.character(final[,1]), 
            ranges = IRanges(as.numeric(final[,2]), as.numeric(final[,3])),
            sites = as.character(final[,4]),
            strand = as.character(final[,5]),
            isCluster = as.logical(final[,6]),
            status = as.character(final[,7]))

            names(final) <- NULL
    }else{
        message("No cluster found")
        return(NULL)
    }

    ##if verbose is FALSE, get only clusters with "PASS"
    if(verbose == FALSE) {
        final = subset(final, final$status == "PASS")
    }

    ##if TRUE, get everything
    if(verbose == TRUE) {
        final
    }

    end.time = Sys.time()
    print(end.time - start.time)
    return(final)
}


load_data <- function(all_files, c) {
    df = NULL
    for(i in seq_along(all_files)){
        if(grepl("\\.bed$", all_files[i])){
            b = import(all_files[i])
        }
        else if(grepl("\\.vcf$|\\.vcf.gz$", all_files[i])){
            b = rowRanges(readVcf(all_files[i]))
        }
        else{
            stop("only .bed, .vcf and .vcf.gz extensions are accepted")
        }
        b$site = names(c[i])
        df = rbind(df, as.data.frame(b)[, c("seqnames", "start", "end",
                                            "strand", "site")])
    }
    df$strand[df$strand == "*"] <- "+"
    return(df)
}


cluster_sites<-function(gr, w, c, overlap, n, res, s, greedy, order,
                        sitesToExclude, sites_orientation, intra_distance){

    start_site <- start(gr)
    end_site <- end(gr)
    end_check <- 0
    site <- gr$site
    exclusion_ls <- sitesToExclude
    strand <- as.vector(strand(gr))
    sites_orientation_input <- sites_orientation
    intraDistance <- intra_distance

    isCluster <- FALSE

    if (greedy == FALSE) {
        upper_boundary = length(gr) - n + 1
    } else {
        upper_boundary = length(gr)
    }
    for (i in seq_len(upper_boundary)) {
        ## check for overlap. basically, the first new
        ## site should satisfy the overlap condition
        if (greedy == TRUE) {
            if (end_site[i] > start_site[i] + w) {
                next
            }
        }
        ## 
        if (isCluster) {
            if ((start_site[i] - end) < overlap) {
                next
            } else {
                isCluster = FALSE
            }
        }
        if (greedy == TRUE) {
            end = end(gr)[end(gr) <= start_site[i] + w & end(gr)>= end_site[i]]
            end = end[length(end)]
            if(end_check == 0){
                end_check = end
                iEnd = which(end(gr) == end)[1]
                ls <- site[i:iEnd]
                so <- strand[i:iEnd]
                sc <- start_site[i:iEnd]
                ec <- end_site[i:iEnd]
            }
            else if(end_check == end & is.null(exclusion_ls)){
                next
            }
            else {
                end_check = end
                iEnd = which(end(gr) == end)[1]
                ls <- site[i:iEnd]
                so <- strand[i:iEnd]
                sc <- start_site[i:iEnd]
                ec <- end_site[i:iEnd]
            }
        }
        else {
            end = max((end_site[i:(i + n - 1)]))
        }
        ## putative cluster is bigger than window
        if (greedy == FALSE) {
            if ((end - start_site[i]) > w) {
                next
            } else {
                ls <- site[i:(i + n - 1)]
                so <- strand[i:(i + n - 1)]
                sc <- start_site[i:(i + n -1)]
                ec <- end_site[i:(i + n -1)]
            }
        }

        ## checking if required sites are present
        ans <- testCombn(ls, c, order, exclusion_ls,so,sites_orientation_input, 
                sc, ec, intraDistance)
        isCluster = ans$logical
        status = ans$status

        ## if we get combnFail, skip
        if(status == "combnFail"){
            next
        }

        res[i, "seqnames"] <- as.character(seqnames(gr)[i])
        res[i, "start"] <- start_site[i]
        res[i, "end"] <- end
        res[i, "strand"] <- s
        res[i, "status"] <- ans$status
        if (greedy == TRUE) {
            res[i, "site"] <- paste(site[i:iEnd], collapse = ",")
        } else{
            res[i, "site"] <- paste(site[i:(i + n - 1)], collapse = ",")
        }
        res[i, "isCluster"] <- isCluster
    }
    
    res <- res[complete.cases(res), ]

    #Case with one entry, Converting vector to matrix prior to return
    if (is(res, "character")){
        res.names <- names(res)
        res = matrix(res, 1, length(res))
        colnames(res) <- res.names
    }
    return(res)
}

testCombn <- function(ls, c, order, sitesToExclude, so, s_orientation_input, 
                    start_site, end_site, intraDistance) {
    ans <- list()
    temp <- NULL
    check <- FALSE
    intra_distance_check <- FALSE


    for (key in names(c)) {
        if (sum((ls == key)) < c[key]) {
            ans$logical = FALSE
            ans$status = "combnFail"
            return(ans)
        }
    }

    if(is.null(order)){ ## order doesn't matter
        ans$logical = TRUE
        ans$status = "PASS"
    }else{
        if(grepl(paste(order,collapse=";"),paste(ls,collapse=";")) == TRUE 
            & is.null(s_orientation_input))
        {
            ans$logical = TRUE
            ans$status = "PASS"
        }
        else if(grepl(paste(order,collapse=";"),paste(ls,collapse=";")) == TRUE 
            & !(is.null(s_orientation_input))){
            for(i in 1:length(ls)){
                temp <- as.vector(na.omit(ls[i:(i+length(order)-1)]))
                if((length(temp) == length(order)) && (all(temp == order)) && 
                    (all(as.vector(na.omit(so[i:(i+length(order)-1)] 
                        == s_orientation_input))))){
                    check <- TRUE
            }}
            if(check){
            ans$logical = TRUE
            ans$status = "PASS"
        }else{
            ans$logical = FALSE
            ans$status = "SitesOrientation"
            }
        }
        else{
            ans$logical = FALSE
            ans$status = "orderFail"
        }
    }
    
    
    if(!(is.null(sitesToExclude))){
        for (exc_site in sitesToExclude){
            if(length(grep(exc_site, ls, value = TRUE))>0){
            ans$logical = FALSE
            ans$status = "ExcludedSites"
        }
        }
    }

    if(intraDistance != 0){
        if(intraDistance > 0){
        #if intraDistance is positive, it means sites should have min distance and above
            for(i in 1:(length(start_site)-1)){
                if((cbind(start_site, end_site)[(i+1),1]) - (cbind(start_site, end_site)[(i),2]) <= intraDistance){
                    intra_distance_check <- TRUE
                    break
                }
            }
        }

        #if intraDistance is negative, it means sites should have max distance and below
        else if(intraDistance < 0){
            for(i in 1:(length(start_site)-1)){
                if((cbind(start_site, end_site)[(i+1),1]) - (cbind(start_site, end_site)[(i),2]) >= abs(intraDistance)){
                    intra_distance_check <- TRUE
                    break
                }
            }
        }
        
        if(intra_distance_check){
            ans$logical = FALSE
            if(ans$status == "PASS"){
            ans$status = "IntraDist"
            }
            else
            {
                ans$status = paste(ans$status, "IntraDist", sep=",")
            }
        }
    }

    
    return(ans)
}

##' @title fcScan

##' @export getCluster

globalVariables(c("parOut","s","site","strand"))

getCluster <- function(x, w, c, overlap = 0, greedy = FALSE, seqnames = NULL,
                    s = "*" , order = NULL, verbose = FALSE) {

    start.time = Sys.time()

    ##check arguments
    if (missing(x)) {
        stop("A data frame, Granges object or input files are needed")
    }

    sitesToExclude <- NULL

    ## Checking input format
    if (is(x, "data.frame") || is(x, "GRanges")) {
        ## print("data.frame input")
        if (is.null(names(c))) {
            stop("When input is a data frame or GRanges object, site names
in condition must be explicitly defined")
        }
        

    } else if (is(x, "GRanges")){
        ## print("GRanges input")
        x <- as.data.frame(x)
    } else {
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
    }
    
    if(length(which(duplicated(names(c))))!=0){
        stop("Names of conditions must be unique")
    }

    if(length(which(c<0))!=0 | !(is.numeric(c))){
        stop("Only positive integers allowed")
    }

    if(!(is.numeric(overlap))){
        stop("Only integers allowed")
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


        # if(greedy == FALSE){
        #     ## order has more sites than condition
        #     for(i in seq_along(unique(order))){
        #         if(sum(unique(order)[i] == order) >
        #             c[which(names(c) == unique(order)[i])]){
        #             stop("Greedy is FALSE and order is larger than condition")
        #         }
        #     }
        # }

    ##check number of sites if equal in condition and order
    if(greedy == FALSE){
        if(!all(sort(c) == summary(as.factor(order)))){
            stop("Greedy is FALSE and order is larger than condition")
            }
        }
    }

    n = sum(c)

    ##getting sites found on the required seqnamesom
    if (length(seqnames) > 0) {
        x = x[(x$seqnames %in% seqnames),]
    }

    ##getting sites found on the required strand
    if (s != "*") {
        x = subset(x, strand %in% s)
    }

    ##check if the sites given in condition c are found in the data
    if( !all(names(c) %in% x$site) ) {
        message("Sites in condition do not match sites in data")
        return(NULL)
    }

    cat (nrow(x), " entries loaded", "\n")

    ##need to subset to keep only the sites required by the user,
    ##if the user wants sites "a","b" and input has "a", "b" and "c",
    ##we subset to get sites "a", "b"
    
    if(!(is.null(sitesToExclude))){
        x
    }
    
    else{
    x = subset(x, site %in% names(c))
    }
    
    unique_seqnames = unique(x$seqnames)

    ## creating an array to fill results
    res = array(data = NA, dim = c(nrow(x), ncol = 7))
    colnames(res) = c("seqnames", "start", "end", "site", "strand",
                        "isCluster", "status")

    final = data.frame()
    
    ## looping over chromosomes 
    for(seq in seq_along(unique_seqnames)){

        df = subset(x, seqnames == unique_seqnames[seq])
        df = df[order(df$start), ]

        if (nrow(df) >= n) {
            result = cluster_sites(df, w, c, overlap, n,
                                    res, s, greedy, order, sitesToExclude)
            final = rbind(final, result)
        }
    }

    row.names(final) <- NULL

    ## if we have only one row the class will be a character
    if (is.character(final)) {
        final = t(final)
    }

    ##check verbose input argument
    if( !(verbose %in% c("TRUE", "FALSE"))) {
        stop("verbose should be TRUE or FALSE")
    }

    ##if verbose is FALSE, get only clusters with "PASS"
    if(verbose == FALSE) {
        final = subset(final, final$status == "PASS")
    }

    ##if TRUE, get everything
    if(verbose == TRUE) {
        final
    }
    
    if (nrow(final) != 0) {
        final = data.frame(seqnames = as.character(final[,"seqnames"]),
                            start = as.numeric(levels(final[,"start"]))
                            [final[,"start"]],
                            end = as.numeric(levels(final[,"end"]))
                            [final[,"end"]],
                            strand = as.character(final[,"strand"]),
                            sites = as.character(final[,"site"]),
                            isCluster = as.logical(final[,"isCluster"]),
                            status = as.character(final[,"status"]),
                            stringsAsFactors = FALSE
                            )

        final <- makeGRangesFromDataFrame(final, keep.extra.columns = TRUE,
                                        starts.in.df.are.0based = TRUE)
    }else{
        message("No cluster found")
        return(NULL)
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
        b$site = names(c[i])
        df = rbind(df, as.data.frame(b)[, c("seqnames", "start", "end",
                                            "strand", "site")])
    }
    df$strand[df$strand == "*"] <- "+"
    return(df)
}


cluster_sites<-function(df, w, c, overlap, n, res, s, greedy, order,
                        sitesToExclude){
    start_site <- df$start
    end_site <- df$end
    site <- df$site
    exclusion_ls <- sitesToExclude

    isCluster <- FALSE

    if (greedy == FALSE) {
        upper_boundary = nrow(df) - n + 1
    } else {
        upper_boundary = nrow(df)
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
            end = df$end[df$end <= start_site[i] + w & df$end >= end_site[i]]
            end = end[length(end)]
            iEnd = which(df$end == end)[1]
            ls <- site[i:iEnd]
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
            }
        }

        ## checking if required sites are present
        ans <- testCombn(ls, c, order, exclusion_ls)
        isCluster = ans$logical
        status = ans$status

        ## if we get combnFail, skip
        if(status == "combnFail"){
            next
        }

        res[i, "seqnames"] <- as.character(df$seqnames[i])
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

testCombn <- function(ls, c, order, sitesToExclude) {
    ans <- list()

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
    } else{
        if(grepl(paste(order,collapse=";"),paste(ls,collapse=";")) == TRUE)
        {
            ans$logical = TRUE
            ans$status = "PASS"
        }else{
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
    
    return(ans)
}

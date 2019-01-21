##' @title fcScan

##' @export getCluster

globalVariables(c("parOut","s","site","strand"))

getCluster <-function(x, w, c, overlap=0, greedy=TRUE, seqnames=NULL,
                      s="." , order=NULL, n_cores=2) {

    start.time = Sys.time()

    ##check arguments
    if (missing(x)) {
        stop("A data frame, Granges object or input files are needed")
    }

    ## Checking input format
    if (is(x, "data.frame")) {
        if (is.null(names(c))) {
            stop("When input is a data frame, names of sites in condition must 
be explicitly defined")
        }
    } else if (is(x, "GRanges")){
        x <- as.data.frame(x)
    } else {
        if (length(c) != length(x))
            stop("Condition and files should be of same length")
        
        ##assigning names to condition if null
        if (is.null(names(c))) {
            names(c) = seq(from = 1,to = length(c),by = 1)
        }
        
        x = load_data(all_files = x, c = c)
    }
    
    
    ## Test if window is a positive integer
    if (w < 0 | w%%1!=0)
        stop("Window size should be a positive integer")
    
    if (!(s %in% c("+", "-", ".")))
        stop('Strand needs to be a valid option. Accepted options are "+", "-"
             and "."')

    ##if not greedy and order is given
    if(is.null(order) == FALSE){
        ##checking for consistency of name in order and condition
        for(i in seq_along(order)){
            if(!order[i] %in% names(c)){
                stop(paste("Site names in order and condition do not match: ", order[i], sep=""))
            }
        }
        if(greedy == FALSE){
            ## order has more sites than condition
            for(i in seq_along(unique(order))){
                if(sum(unique(order)[i] == order) > c[which(names(c) == unique(order)[i])]){
                    stop("Greedy is FALSE and order is larger than condition")
                }
            }
        }
    }

    n = sum(c)
    res = array(data = NA, dim = c(nrow(x), ncol = 8))
    colnames(res) = c("seqnames", "start", "end", "size", "site", "strand",
                      "isCluster", "status")
    ##getting sites found on the required seqnamesom
    if (length(seqnames) > 0) {
        x = x[(x$seqnames %in% seqnames),]
    }

    ##getting sites found on the required strand
    if (s != ".") {
        x = subset(x, strand %in% s)
    }

    ##check if the sites given in condition c are found in the data
    if( !all(names(c) %in% x$site) ) {
        message("Sites in condition do not match sites in data")
        return(NULL)
    }

    ##need to subset to keep only the sites required by the user,
    ##if the user wants sites "a","b" and input has "a", "b" and "c", we subset to get sites "a", "b"
    x = subset(x, site %in% names(c))

    unique_seqnames = unique(x$seqnames)

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    final = foreach(parOut = seq_along(unique_seqnames) ,
                    .export = c("testCombn","cluster_sites"),
                    .combine = rbind) %dopar% {
                        df = subset(x, seqnames == unique_seqnames[parOut])
                        df = df[order(df$start), ]
                        if (nrow(df) >= n) {
                            result = cluster_sites(df, w, c, overlap, n,
                                                   res, s, greedy, order)
                        }
                    }
    row.names(final) <- NULL

    ## if we have only one row the class will be a character
    if (class(final) == "character") {
        final = t(final)
    }
    
    if (nrow(final) != 0) {
        final = data.frame(seqnames = as.character(final[,"seqnames"]),
                           start = as.numeric(final[,"start"]),
                           end = as.numeric(final[,"end"]),
                           size = as.numeric(final[,"size"]),
                           isCluster = as.logical(final[,"isCluster"]),
                           id = paste("c", seq.int(nrow(final)), sep = ""),
                           status = as.character(final[,"status"]),
                           sites = as.character(final[,"site"]),
                           score = 1,
                           stringsAsFactors = FALSE
                           )
    }
    
    final <- makeGRangesFromDataFrame(final, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
    
    stopCluster(cl)
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
        df = rbind(df, as.data.frame(b)[, c("seqnames","start","end","strand","site")])
    }
    df$strand[df$strand == "*"] <- "+"
    return(df)
}


cluster_sites <-function(df, w, c, overlap, n, res, s, greedy, order) {
    start <- df$start
    e <- df$end
    site <- df$site

    isCluster <- FALSE

    if (greedy == FALSE) {
        upper_boundary = nrow(df) - n + 1
    } else {
        upper_boundary = nrow(df)
    }

    ## non greedy
    for (i in seq_along(1:upper_boundary)) {
        ## check for overlap. basically, the first new
        ##site should satisfythe overlap condition
        if (greedy == TRUE) {
            if (e[i] > start[i] + w) {
                next
            }
        }
        if (isCluster) {
            if ((start[i] - end) < overlap) {
                next
            } else {
                isCluster = FALSE
            }
        }
        if (greedy == TRUE) {
            end = df$end[df$end <= start[i] + w & df$end >= e[i]]
            end = end[length(end)]
            iEnd = which(df$end == end)[1]
            ls <- site[i:iEnd]
        }
        else {
            end = max((e[i:(i + n - 1)]))
        }
        ## putative cluster is bigger than window
        if (greedy == FALSE) {
            if ((end - start[i]) > w) {
                next
            } else {
                ls <- site[i:(i + n - 1)]
            }
        }

        ## checking if required sites are present
        ans <- testCombn(ls, c, order)
        isCluster = ans$logical

        res[i, "seqnames"] <- as.character(df$seqnames[i])
        res[i, "start"] <- start[i]
        res[i, "end"] <- end
        res[i, "size"] <- end - start[i]
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
    return(res)
}


testCombn <- function(ls, c, order) {
    ans <- list()

    for (key in names(c)) {
        if (sum((ls == key)) < c[key]) {
            ans$logical = FALSE
            ans$status = "combnFail"
            return(ans)
        }
    }

    if(is.null(order)){ ##order doesnt matter
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
    return(ans)
}




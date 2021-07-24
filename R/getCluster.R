##' @title fcScan

##' @export getCluster
globalVariables(c("chr"))
getCluster <- function(x, w, c, overlap = 0, greedy = FALSE, seqnames = NULL,
                    s = "*" , order = NULL, site_orientation = NULL, 
                    site_overlap = 0, cluster_by = "startsEnds", 
                    allow_clusters_overlap = FALSE, 
                    include_partial_sites = FALSE, 
                    partial_overlap_percentage= NULL, 
                    threads = 1, verbose = FALSE) {

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
        stop("Only positive integers are allowed")
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

    ##check if number of sites is equal betweeb condition and order
    if(greedy == FALSE){
        count_elements <- c(count(order)$freq)
        names(count_elements) <- count(order)$x
        if(!(all(sort(count_elements) == sort(c)))){
            stop("Greedy is FALSE and order is smaller/larger than condition")
            }
        }
    }
    
    #site orientation works only if order is provided
    if(!(is.null(site_orientation)) & is.null(order)){
        stop("site_orientation cannot be used unless order is specified")
    }

    #site orientation should be positive or negative
    if(!(all(site_orientation %in% c("+","-")))){
        stop('Site orientation should be "+" or "-"')
    }

    #site orientation input length should have same order length
    if(length(order) != length(site_orientation) & 
        !(is.null(site_orientation))){
        stop("When specified, sites orientation must be defined to
        all sites following 'order' option")
    }

    ##site_overlap should be either positive, negative or zero
    if(!(is.numeric(site_overlap) || site_overlap == 0)){
        stop("Only integers are allowed for distance between sites")
    }

    ##check cluster_by input
    if(!(cluster_by %in% c("startsEnds","endsStarts","starts",
            "ends","middles")) || length(cluster_by)>1 || is.null(cluster_by)){
        stop("Invalid cluster_by option")
    }

    ##check allow_clusters_overlap usage
    #if(is.null(cluster_by) && allow_clusters_overlap == TRUE){
    #    warning("Use 'allow_clusters_overlap' when 'cluster_by' is assigned")
    #}

    ##check include_partial_sites usage (1)
    #if(is.null(cluster_by) && include_partial_sites == TRUE){
    #    warning("Use 'include_partial_sites' when 'cluster_by' is assigned")
    #}
    
    ##check include_partial_sites usage (2)
    if(!(cluster_by %in% c("startsEnds","ends")) 
        && include_partial_sites == TRUE){
        warning("Use 'include_partial_sites' only with 'startsEnds' or 'ends'")
    }

    ##check partial_overlap_percentage (1)
    if(!is.null(partial_overlap_percentage) && 
        !is.numeric(partial_overlap_percentage)){
            stop("Invalid partial_overlap_percentage input")
        }
            
    ##check partial_overlap_percentage (2)
    if(is.numeric(partial_overlap_percentage) && 
        (length(partial_overlap_percentage)>1 || 
        partial_overlap_percentage > 1 || 
        partial_overlap_percentage < 0)){
            stop("Invalid partial_overlap_percentage input")
        }

    ##check if partial_overlap_percentage equals zero
    if(!is.null(partial_overlap_percentage) && 
        partial_overlap_percentage == 0){
        partial_overlap_percentage = NULL
    }

    ##check if partial_overlap_percentage is used without include_partial_sites
    if(!(is.null(partial_overlap_percentage)) && !include_partial_sites){
        warning("partial_overlap_percentage ignored - include_partial_sites=F")
    }

    ##check verbose input argument
    if( !(verbose %in% c("TRUE", "FALSE"))) {
        stop("Verbose should be TRUE or FALSE")
    }

    ## n contains the total number of sites desired
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

    if(threads %% 1 !=0 || threads < 0 || threads > detectCores()){
        message("invalid \"threads\" value")
        message("getCluster will use all the available threads instead")
        threads = 0
    }
    if(threads == 0){
        if(detectCores() >= length(unique(as.character(seqnames(x))))){
            threads = length(unique(as.character(seqnames(x))))
        }
        else if(detectCores() < length(unique(as.character(seqnames(x))))){
            threads = detectCores()
        }
    } 
    if(threads > length(unique(as.character(seqnames(x))))){
        threads = length(unique(as.character(seqnames(x))))
        cat("Some threads are not going to be used\n")
    }
    
    unique_seqnames = unique((seqnames(x)))

    ## creating an array to fill results
    res = array(data = NA, dim = c(length(x), ncol = 7))
    colnames(res) = c("seqnames", "start", "end", "site", "strand",
                        "isCluster", "status")


    ## looping over chromosomes and call cluster_sites
    cat("Running getCluster using",threads,"threads\n")
    registerDoParallel(threads)
    if(length(x) >= n && greedy == FALSE){
        final <- foreach(chr=unique(as.character(seqnames(x))),
        .combine = rbind) %dopar% {
            gr = subset(x, seqnames == chr)
            gr = sort(gr, ignore.strand = TRUE)
            cluster_by_greedy_false(gr, w, c, overlap, n,
                    res, s, greedy, order, sitesToExclude, 
                    site_orientation, site_overlap, allow_clusters_overlap,
                    cluster_by, include_partial_sites, 
                    partial_overlap_percentage)
        }
    } else if(length(x) >= n && greedy == TRUE){
        final <- foreach(chr=unique(as.character(seqnames(x))),
        .combine = rbind) %dopar% {
            gr = subset(x, seqnames == chr)
            gr = sort(gr, ignore.strand = TRUE)
            cluster_by_greedy_true(gr, w, c, overlap, n,
                    res, s, greedy, order, sitesToExclude, 
                    site_orientation, site_overlap, allow_clusters_overlap,
                    cluster_by, include_partial_sites, 
                    partial_overlap_percentage)
        }
    }

    
    
    # for(seq in seq_along(unique_seqnames)){

    #     gr = subset(x, seqnames == unique_seqnames[seq])
    #     gr = sort(gr, ignore.strand = TRUE)

    #     if (length(gr) >= n) {
    #         #cluster_by is NULL - default algorithm is used
    #         # if(is.null(cluster_by)){
    #         # result = cluster_sites(gr, w, c, overlap, n,
    #         #         res, s, greedy, order, sitesToExclude, 
    #         #             site_orientation, site_overlap)
    #         # final = rbind(final, result)
    #         # }

    #         #valid cluster_by option and greedy = FALSE
    #         if(greedy == FALSE){
    #             result = cluster_by_greedy_false(gr, w, c, overlap, n,
    #                 res, s, greedy, order, sitesToExclude, 
    #                 site_orientation, site_overlap, allow_clusters_overlap,
    #                 cluster_by, include_partial_sites, 
    #                 partial_overlap_percentage)
    #         final = rbind(final, result)
    #         }

    #         #valid cluster_by option and greedy = TRUE
    #         else if(greedy == TRUE){
    #             result = cluster_by_greedy_true(gr, w, c, overlap, n,
    #                 res, s, greedy, order, sitesToExclude, 
    #                 site_orientation, site_overlap, allow_clusters_overlap,
    #                 cluster_by, include_partial_sites, 
    #                 partial_overlap_percentage)
    #         final = rbind(final, result)
    #         }
    #     }
    # }

    #Get final result
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
    #Case when input is BED and/or VCF files
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

## n contains the total number of sites desired
## s contains strand
## res is array for temporary results

cluster_by_greedy_false <- function(gr, w, c, overlap, n, res, s, greedy, order,
                        sitesToExclude, site_orientation, site_overlap,
                        allow_clusters_overlap, cluster_by, 
                        include_partial_sites, partial_overlap_percentage){

    #Greedy False - different clustering options
    #Argument: cluster_by
    #Description: clustering sites using different options.
    #cluster_by options: startEnds, endsStarts, starts, ends, middles

    #Argument: include_partial_sites
    #Description: this argument allows to include partially overlapping sites
    # with window size at the end of the cluster. If set to 'TRUE', partially
    # overlapping sites are included, however, if set to 'FALSE', overlapping
    # sites at the end of the cluster are not included. This option is used
    # only when 'cluster_by' is either 'startsEnds' or 'ends'
    # allow_clusters_overlap gives the user the ability to get overlapping
    # clusters or not.

    #Sub-clusters are omitted 
    #(clusters having different start coordinate but same end coordinate)

    start_site <- start(gr)
    end_site <- end(gr)
    end_check <- 0
    site <- gr$site
    exclusion_ls <- sitesToExclude
    strand <- as.vector(strand(gr))
    site_orientation_input <- site_orientation
    site_overlap_input <- site_overlap

    isCluster <- FALSE
    
    #looping index - esseential for jumping clusters or sites
    #it depends on user input: allow_clusters_overlap T or F

    for(i in seq_len(length(gr)-n+1)){
        #Check overlap - Zero overlap returns all clusters
        #either overlapping or not
        if (isCluster) {
            if ((start_site[i] - end) < overlap & overlap !=0) {
                next
                ## Allow or deny overlapping clusters
            } else if((start_site[i] - end) < overlap & 
                        allow_clusters_overlap == FALSE){
                next
            } else {
                isCluster = FALSE
            }
        }

        #cluster_by: startEnds
        if(cluster_by == "startsEnds"){
        if(!include_partial_sites){
            if(end_site[(i+n-1)] - start_site[i] <= w){
                ls <- site[i:((i+n)-1)]
                sc <- start_site[i:((i+n)-1)]
                ec <- end_site[i:((i+n)-1)]
                so <- strand[i:((i+n)-1)]
                end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

                end <- ec[length(ec)]
            } else {
                next
            }
        }

        if(include_partial_sites){
            if (is.null(partial_overlap_percentage)) {
                if((end_site[(i+n-1)] - start_site[i] <= w) || 
                    ((end_site[(i+n-1)] > w + start_site[i]) && 
                    (start_site[i] + w > start_site[(i+n-1)]) && 
                    (start_site[i] + w < end_site[(i+n-1)]))){
                
                    ls <- site[i:((i+n)-1)]
                    sc <- start_site[i:((i+n)-1)]
                    ec <- end_site[i:((i+n)-1)]
                    so <- strand[i:((i+n)-1)]
                    end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

                    end <- ec[length(ec)]
                }
                else {
                    next
                }
            } else { # else if partial_overlap_percentage is set
                if((end_site[(i+n-1)] - start_site[i] <= w) || 
                    ((end_site[(i+n-1)] > w + start_site[i]) && 
                    (start_site[i] + w > start_site[(i+n-1)]) && 
                    (start_site[i] + w < end_site[(i+n-1)])) && 
                    (((end_site[(i+n-1)] - start_site[(i+n-1)])-
                    (end_site[(i+n-1)]-(w+start_site[i])))/(end_site[(i+n-1)]
                    - start_site[(i+n-1)])>=partial_overlap_percentage)){
                
                    ls <- site[i:((i+n)-1)]
                    sc <- start_site[i:((i+n)-1)]
                    ec <- end_site[i:((i+n)-1)]
                    so <- strand[i:((i+n)-1)]
                    end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

                    end <- ec[length(ec)]
                }
                else {
                    next
                }
            }
        }

        #cluster_by: endsStarts
    } else if(cluster_by == "endsStarts"){
        if(start_site[(i+n-1)] - end_site[i] <= w){
            ls <- site[i:((i+n)-1)]
            sc <- start_site[i:((i+n)-1)]
            ec <- end_site[i:((i+n)-1)]
            so <- strand[i:((i+n)-1)]
            end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

            end <- ec[length(ec)]
        } else {
            next
        }

        #cluster_by: middles
    } else if(cluster_by == "middles"){
        if((start_site[(i+n-1)] + floor(width(gr)[(i+n-1)]/2)) - 
            (start_site[i] + floor(width(gr)[i]/2)) 
                <= w){
            ls <- site[i:((i+n)-1)]
            sc <- start_site[i:((i+n)-1)]
            ec <- end_site[i:((i+n)-1)]
            so <- strand[i:((i+n)-1)]
            end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

            end <- ec[length(ec)]
        } else {
            next
        }

        #cluster_by: starts
    } else if(cluster_by == "starts"){
        if(start_site[(i+n-1)] - start_site[i] <= w){
            ls <- site[i:((i+n)-1)]
            sc <- start_site[i:((i+n)-1)]
            ec <- end_site[i:((i+n)-1)]
            so <- strand[i:((i+n)-1)]
            end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

            end <- ec[length(ec)]
        } else {
        next
        }

        #cluster_by: ends
    } else if(cluster_by == "ends"){
        if(!include_partial_sites && is.null(partial_overlap_percentage)){
            if(end_site[(i+n-1)] - end_site[i] <= w){
                ls <- site[i:((i+n)-1)]
                sc <- start_site[i:((i+n)-1)]
                ec <- end_site[i:((i+n)-1)]
                so <- strand[i:((i+n)-1)]
                end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

                end <- ec[length(ec)]
            } else {
                next
            }
        }

        if(include_partial_sites && is.null(partial_overlap_percentage)){
            if((end_site[(i+n-1)] - end_site[i] <= w) || 
                (end_site[(i+n-1)] > w + end_site[i]) && 
                (end_site[i] + w > start_site[(i+n-1)]) && 
                (end_site[i] + w < end_site[(i+n-1)])){

                    ls <- site[i:((i+n)-1)]
                    sc <- start_site[i:((i+n)-1)]
                    ec <- end_site[i:((i+n)-1)]
                    so <- strand[i:((i+n)-1)]
                    end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

                    end <- ec[length(ec)]
            } else {
            next
            }
        }

        if(include_partial_sites && !(is.null(partial_overlap_percentage))){
            if((end_site[(i+n-1)] - end_site[i] <= w) || 
                (end_site[(i+n-1)] > w + end_site[i]) && 
                (end_site[i] + w > start_site[(i+n-1)]) && 
                (end_site[i] + w < end_site[(i+n-1)]) &&
                (((end_site[(i+n-1)] - start_site[(i+n-1)])-
                (end_site[(i+n-1)]-(w+end_site[i])))/(end_site[(i+n-1)]
                - start_site[(i+n-1)]))>=partial_overlap_percentage){
            
                    ls <- site[i:((i+n)-1)]
                    sc <- start_site[i:((i+n)-1)]
                    ec <- end_site[i:((i+n)-1)]
                    so <- strand[i:((i+n)-1)]
                    end_check <- end_site[i:((i+n)-1)][((i+n)-1)]

                    end <- ec[length(ec)]
            } else {
            next
            }
        }
    }

    #test the cluster of sites for different parameters
    ans <- testCombn(ls, c, order, exclusion_ls,so,site_orientation_input, 
                sc, ec, site_overlap_input)
        isCluster = ans$logical
        status = ans$status

        ## if we get combnFail, skip
        if(status == "combnFail"){
            next
        }

        res[i, "seqnames"] <- as.character(seqnames(gr)[i])
        res[i, "start"] <- start_site[i]
        res[i, "end"] <- max(end_site[i:((i+n)-1)])
        res[i, "strand"] <- s
        res[i, "status"] <- status
        res[i, "site"] <- paste(ls, collapse = ",")
        res[i, "isCluster"] <- isCluster
    }
    
    #return resulting cluster - either true or false cluster
    res <- res[complete.cases(res), ]

    #Case with one entry, Converting vector to matrix prior to return
    if (is(res, "character")){
        res.names <- names(res)
        res = matrix(res, 1, length(res))
        colnames(res) <- res.names
    }
    return(res)
}

cluster_by_greedy_true <- function(gr, w, c, overlap, n, res, s, greedy, order,
                        sitesToExclude, site_orientation, site_overlap, 
                        allow_clusters_overlap, cluster_by, 
                        include_partial_sites, partial_overlap_percentage){

    #Greedy True - different clustering options
    #Argument: cluster_by
    #Description: clustering sites using different options.
    #cluster_by options: startEnds, endsStarts, starts, ends, middles

    #Argument: include_partial_sites
    #Description: this argument allows to include partially overlapping sites
    # with window size at the end of the cluster. If set to 'TRUE', partially 
    # overlapping sites are included, however, if set to 'FALSE', overlapping
    # sites at the end of the cluster are not included. This option is used 
    # only when 'cluster_by' is either 'startsEnds' or 'ends'
    # allow_clusters_overlap gives the user the ability to get overlapping 
    # clusters or not.
    
    #Sub-clusters are omitted 
    #(clusters having different start coordinate but same end coordinate)

    start_site <- start(gr)
    end_site <- end(gr)
    end_coordinate_check <- 0
    site <- gr$site
    exclusion_ls <- sitesToExclude
    strand <- as.vector(strand(gr))
    site_orientation_input <- site_orientation
    site_overlap_input <- site_overlap

    isCluster <- FALSE

    for(i in seq_len(length(gr)-1)){
        #Check overlap - Zero overlap returns all clusters
        #either overlapping or not
        if (isCluster) {
            if ((start_site[i] - end) < overlap & overlap !=0) {
                next
                ## Allow or deny overlapping clusters
            } else if((start_site[i] - end) < overlap & 
                        allow_clusters_overlap == FALSE){
                next
            } else {
                isCluster = FALSE
            }
        }

        #A vector that hold indices of found sites within window size    
        indices <- NULL

        #cluster_by: startEnds
        if(cluster_by == "startsEnds"){
            indices <- which(end_site <= start_site[i] + w)[
                i:(length(end_site <= start_site[i] + w))]
            indices <- indices[!is.na(indices)]

            if(include_partial_sites && is.null(partial_overlap_percentage)){
                last_site <- indices[length(indices)]
                while(last_site < length(gr)){
                    if(end_site[(last_site+1)] >= start_site[i] + w & 
                    start_site[i] + w >= start_site[(last_site+1)] & 
                    start_site[i] + w <= end_site[(last_site+1)]){
                    indices <- c(indices, (last_site+1))
                    indices <- indices[!is.na(indices)]
                    } else {break}
                    last_site = last_site + 1
                }
            }

            if(include_partial_sites && !(is.null(partial_overlap_percentage))){
                last_site <- indices[length(indices)]
                while(last_site < length(gr)){
                    if(end_site[(last_site+1)] >= start_site[i] + w & 
                    start_site[i] + w >= start_site[(last_site+1)] & 
                    start_site[i] + w <= end_site[(last_site+1)] & 
                    (((end_site[(last_site+1)] - start_site[(last_site+1)])-
                    (end_site[(last_site+1)]-(w+start_site[i])))/
                    (end_site[(last_site+1)]- start_site[(last_site+1)]))
                        >=partial_overlap_percentage){
                    indices <- c(indices, (last_site+1))
                    indices <- indices[!is.na(indices)]
                    } else {break}
                    last_site = last_site + 1
                }
            }
        }

        #cluster_by: endsStarts
        else if(cluster_by == "endsStarts"){
            indices <- which(start_site <= end_site[i] + w)[
                i:(length(start_site <= end_site[i] + w))]
            indices <- indices[!is.na(indices)]
        }

        #cluster_by: middles        
        else if(cluster_by == "middles"){
            indices <- which(start_site + floor(width(gr)/2) <= 
                (start_site[i] + floor(width(gr)[i]/2)) + w)[
                    i:length(start_site + floor(width(gr)/2) <= 
                        (start_site[i] + floor(width(gr)[i]/2)) + w)]
            indices <- indices[!is.na(indices)]
        }

        #cluster_by: starts
        else if(cluster_by == "starts"){
            indices <- which(start_site <= start_site[i] + w)[
                i:(length(start_site <= start_site[i] + w))]
            indices <- indices[!is.na(indices)]
        }

        #cluster_by: ends
        else if(cluster_by == "ends"){
            indices <- which(end_site <= end_site[i] + w)[
                i:(length(end_site <= end_site[i] + w))]
            indices <- indices[!is.na(indices)]
            
            if(include_partial_sites && is.null(partial_overlap_percentage)){
                last_site <- indices[length(indices)]
                while(last_site < length(gr)){
                if((end_site[last_site+1] - end_site[i] > w) && 
                    (end_site[i] + w >= start_site[last_site+1]) && 
                    (end_site[i] + w <= end_site[last_site+1])){
                            indices <- c(indices, last_site+1)
                } else {break}
                last_site = last_site + 1
                }
            }

            else if(include_partial_sites && 
            !(is.null(partial_overlap_percentage))){
                last_site <- indices[length(indices)]
                while(last_site < length(gr)){
                if((end_site[last_site+1] - end_site[i] > w) && 
                    (end_site[i] + w >= start_site[last_site+1]) && 
                    (end_site[i] + w <= end_site[last_site+1]) && 
                    (((end_site[(last_site+1)] - start_site[(last_site+1)])-
                    (end_site[(last_site+1)]-(w+end_site[i])))/
                    (end_site[(last_site+1)] - start_site[(last_site+1)]))
                        >=partial_overlap_percentage){
                            indices <- c(indices, last_site+1)
                } else {break}
                last_site = last_site + 1
                }
            }
        }
        
        #get 'end_site' coordinate
        end <- end_site[indices[length(indices)]]

        #Jump cluster if it has less than number of minimum
        #required sites
        if(length(indices) < n){
            next
        } else if(end_coordinate_check == 0) { 
            #First cluster does not require sub-cluster check
            ls <- gr$site[indices]
            so <- as.vector(strand(gr)[indices])
            sc <- start_site[indices]
            ec <- end_site[indices]
            end_coordinate_check = end
        } else if (end_coordinate_check == end && is.null(exclusion_ls)){
            next #omit sub-cluster
        } else if(end_coordinate_check == end && !(is.null(exclusion_ls))){
            ls <- gr$site[indices]
            so <- as.vector(strand(gr)[indices])
            sc <- start_site[indices]
            ec <- end_site[indices]
            end_coordinate_check = end
        } else if(end_coordinate_check != end){
            ls <- gr$site[indices]
            so <- as.vector(strand(gr)[indices])
            sc <- start_site[indices]
            ec <- end_site[indices]
            end_coordinate_check = end
        }

    #test the cluster of sites for different parameters
    ans <- testCombn(ls, c, order, exclusion_ls,so,site_orientation_input, 
                sc, ec, site_overlap_input)
        isCluster = ans$logical
        status = ans$status

        ## if we get combnFail, skip
        if(status == "combnFail"){
            next
        }

        res[i, "seqnames"] <- as.character(seqnames(gr)[i])
        res[i, "start"] <- start_site[i]
        res[i, "end"] <- max(end_site[indices])
        res[i, "strand"] <- s
        res[i, "status"] <- status
        res[i, "site"] <- paste(ls, collapse = ",")
        res[i, "isCluster"] <- isCluster
    }
    
    #return resulting cluster - either true or false cluster
    res <- res[complete.cases(res), ]

    #Case with one entry, Converting vector to matrix prior to return
    if (is(res, "character")){
        res.names <- names(res)
        res = matrix(res, 1, length(res))
        colnames(res) <- res.names
    }
    return(res)
}

testCombn <- function(ls, c, order, sitesToExclude, so, s_orientation, 
                    start_site, end_site, s_overlap) {
    ans <- list()
    temp <- NULL
    site_orientation_check <- FALSE
    site_overlap_check <- FALSE

    #check for minimum number of sites if satisfying condition
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
    }else{ # check for order validity
        if(grepl(paste(order,collapse=";"),paste(ls,collapse=";")) == TRUE 
            & is.null(s_orientation))
        {
            ans$logical = TRUE
            ans$status = "PASS"
        }
        # check for sites orientation validity 
        else if(grepl(paste(order,collapse=";"),paste(ls,collapse=";")) == TRUE 
            & !(is.null(s_orientation))){
            for(i in seq_along(ls)){#use of seq_along as advised by BiocCheck
                temp <- as.vector(na.omit(ls[i:(i+length(order)-1)]))
                if((length(temp) == length(order)) && (all(temp == order)) && 
                    (all(as.vector(na.omit(so[i:(i+length(order)-1)] 
                        == s_orientation))))){
                    site_orientation_check <- TRUE
            }}
            if(site_orientation_check){
            ans$logical = TRUE
            ans$status = "PASS"
        }else{
            ans$logical = FALSE
            ans$status = "siteOrientation"
            }
        }
        else{
            ans$logical = FALSE
            ans$status = "orderFail"
        }
    }
    
    # check for excluded sites - zero sites in condition
    if(!(is.null(sitesToExclude))){
        for (exc_site in sitesToExclude){
            if(length(grep(exc_site, ls, value = TRUE))>0){
            ans$logical = FALSE
            ans$status = "excludedSites"
        }
        }
    }

    # check for distance between sites within cluster
    if(s_overlap != 0){
        if(s_overlap > 0){
        #if s_overlap is positive, 
        #it means sites should have min distance and above
            for(i in seq_len(length(start_site)-1)){
                if((cbind(start_site, end_site)[(i+1),1]) - 
                    (cbind(start_site, end_site)[(i),2]) < s_overlap){
                    site_overlap_check <- TRUE
                    break
                }
            }
        }

        #if s_overlap is negative, 
        #it means sites should have max distance and below
        else if(s_overlap < 0){
            for(i in seq_len(length(start_site)-1)){
                if((cbind(start_site, end_site)[(i+1),1]) - 
                    (cbind(start_site, end_site)[(i),2]) > abs(s_overlap)){
                    site_overlap_check <- TRUE
                    break
                }
            }
        }
        
        if(site_overlap_check){
            if(ans$status == "PASS"){
            ans$status = "siteOverlap"
            ans$logical = FALSE
            }
        }
    }
    return(ans)
}
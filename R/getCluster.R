##' @title fcScan

##' @export getCluster

##   Build and Reload Package:  'Ctrl + Shift + B'
##   Check Package:             'Ctrl + Shift + E'
##   Test Package:              'Ctrl + Shift + T'

globalVariables(c("parOut","s","site","strand"))

getCluster <-function(x, w, c, overlap=0, greedy=TRUE, chr=NULL,
                      s="." , order=NULL, n_cores=2) {

    start.time = Sys.time()

    ##check arguments
    if (missing(x)) {
        stop("Need to specify the files to be used")
    }

    ##given data frame and not correct format of condition
    if (is.data.frame(x)) {
        if (is.null(names(c))) {
            stop("When input is data frame condition must be explicitly
                 defined")
        }
        } else {
            if (length(c) != length(x))
                stop("Condition and files should be of same length")
            ##assigning condition
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
        ##checking same consistency of name in order and condition
        for(i in 1:length(order)){
            if(!(order[i] %in% names(c))){
                stop("Site names in order and condition do not match")
            }
        }
        if(greedy == FALSE){
            if(length(order) > sum(c)) ##order has more sites than condition
                stop("Order has more sites than condition")
        }
    }

    n = sum(c)
    res = array(data = NA, dim = c(nrow(x), ncol = 8))
    colnames(res) = c("chr", "start", "end", "size", "site", "strand",
                      "isCluster", "status")
    ##getting sites found on the required chrom
    if (length(chr) > 0) {
        x = x[(x$chr %in% chr),]
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
    ##if the user wants sites 1,2 and the file has 1,2,3 we subset to get sites 1,2
    x = subset(x, site %in% names(c))

    unique_chr = unique(x$chr)

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    final = foreach(parOut = 1:length(unique_chr) ,
                    .export = c("testCombn","cluster_sites"),
                    .combine = rbind) %dopar% {
                        df = subset(x, chr == unique_chr[parOut])
                        df = df[order(df$start), ]
                        if (nrow(df) >= n) {
                            result = cluster_sites(df, w, c, overlap, n,
                                                   res, s, greedy, order)
                        }
                    }
    row.names(final) <- NULL

    ##if we have only one row the class will be a character
    if (class(final) == "character") {
        final = t(final)
    }
    final = as.data.frame(final , stringsAsFactors = FALSE)
    if (nrow(final) != 0) {
        final$start = as.numeric(final$start)
        final$end = as.numeric(final$end)
        final$size = as.numeric(final$size)
        final$isCluster = as.logical(final$isCluster)
        final$id = paste("c", seq.int(nrow(final)), sep = "")
        final$score = 1
    }
    stopCluster(cl)
    end.time = Sys.time()
    print(end.time - start.time)
    return(final)
}



load_data <- function(all_files, c) {
    bed_files = c()
    vcf_files = c()

    ## getting bed files into bed_files and vcf files into vcf_files
    bed_files = grep("\\.bed$", all_files, value = TRUE)
    vcf_files = grep("\\.vcf$|\\.vcf.gz$", all_files, value = TRUE)

    if (length(bed_files) != 0) {
        data_bed = load_files(files=bed_files, c=c, x=all_files)
    }

    if (length(vcf_files) != 0) {
        data_vcf = load_files(files=vcf_files, c=c, x=all_files)
    }

    if (length(bed_files) != 0 & length(vcf_files) != 0) { #we have bed and vcf files
        data <-  rbind(data_bed, data_vcf)
        rownames(data) <- c()
        return(data)
    }
    else if (length(bed_files) != 0 & length(vcf_files) == 0) {#we have bed files
        rownames(data_bed) <- c()
        return(data_bed)
    }
    else {#we have vcf files
        rownames(data_vcf) <- c()
        return(data_vcf)
    }
}

load_files <- function(files, c, x) {
    site = c()
    print (files)
    if(file_ext(files[1]) == "bed") {
        data_f <- lapply(files, import)
    } else {
        data_f <- lapply(files, readVcf)
    }

    names(data_f) <- files

    ## extract the needed cols: names,start,end,strand,site

    chr = unlist(lapply(data_f, function(x) as.character(seqnames(x))), use.names = FALSE)
    start_sites = unlist(lapply(data_f, function(x) start(x)), use.names = FALSE)
    end_sites = unlist(lapply(data_f, function(x) end(x)), use.names = FALSE)
    strand = unlist(lapply(data_f, function(x) as.character(strand(x))), use.names = FALSE)

    for(i in seq_along(files)){
        index = which(x == files[i])
        site = c(site, rep(names(c[index]), length(data_f[[i]])))
    }

    ## creating the dataframe from granges
    df1 = data.frame("chr"= chr ,"start" = start_sites,
                     "end" = end_sites,"strand"= strand,
                     "site" = site, stringsAsFactors = FALSE)

    if(length(which(df1$strand == "*")) > 0){
        df1[which(df1$strand == "*"),]$strand = "+"
    }
    print("&***********************")
    return(df1)
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
    for (i in 1:upper_boundary) {
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

        res[i, "chr"] <- as.character(df$chr[i])
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




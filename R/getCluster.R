#' @title fcScan

#' @export getCluster

#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

globalVariables(c("parOut","s","site","strand"))

getCluster <-function(x, w, c, overlap=0, greedy=TRUE, chr=NULL,
                      s="." , n_cores=2, order=NULL) {

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
            x = load_data(my_files = x, c = c)
    }

    if (w < 0)
        stop("Window size cannot be a negative number")

    if (!(s %in% c("+", "-", ".")))
        stop('Strand needs to be a valid option. Accepted options are "+", "-"
     and "."')

    ##if not greedy and order is given

        if(is.null(order) == FALSE){
            ##checking same consistency of name in order and condition
            for(i in 1:length(order)){
                if(!(order[i] %in% names(c))){
                    stop("Site names in order and condition not consistent")
                }
            }
            if(greedy == FALSE){
            if(length(order) > sum(c)) ##order has more sites than condition
                stop("Not convenient order and condition")
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

    ##need to subset to keep only the sites required by the user
    x = subset(x, site %in% names(c))
    if (nrow(x) == 0) {
        message("Required sites not found")
        stopCluster(cl)
        return(NULL)
    }

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



load_data <- function(my_files, c) {
    bed_files = c()
    vcf_files = c()
    ##adding files to vcf or bed depending on extension
    bed_files = grep("\\.bed$", my_files, value = TRUE)
    vcf_files = grep("\\.vcf$|\\.vcf.gz$", my_files, value = TRUE)

    if (length(bed_files) != 0) {
        data_bed = load_bed_files(bed_files=bed_files, c=c, x=my_files)
    }

    if (length(vcf_files) != 0) {
        data_vcf = load_vcf_files(vcf_files=vcf_files, c=c, x=my_files)
    }

    if (length(bed_files) != 0 & length(vcf_files) != 0) {
        data <- do.call("rbind", c(data_bed, data_vcf))
        rownames(data) <- c()
        return(data)
    }
    else if (length(bed_files) != 0 & length(vcf_files) == 0) {
        data <- do.call("rbind", (data_bed))
        rownames(data) <- c()
        return(data)
    }
    else {
        data <- do.call("rbind", (data_vcf))
        rownames(data) <- c()
        return(data)
    }
}

load_bed_files <- function(bed_files, c, x) {
    data_bed <- lapply(bed_files, read.table, stringsAsFactors = FALSE)

    names(data_bed) <- bed_files
    for (i in seq_along(data_bed)) {
        if (ncol(data_bed[[i]]) >= 6) {
            ##extrat col(1,2,3,6) has strand col
            data_bed[[i]] = data_bed[[i]][, c(1, 2, 3, 4, 6)]
            colnames(data_bed[[i]]) = c("chr", "start", "end", "name", "strand")
            ##check if strand is "." and replace with +
            if (all(data_bed[[i]]$strand == ".")) {
                data_bed[[i]]$strand = "+"
            }
        } else if (ncol(data_bed[[i]]) == 4 |
                   ncol(data_bed[[i]]) == 5) {
            ##extract col(1,2,3) and add strand col
            data_bed[[i]] = data_bed[[i]][, c(1, 2, 3, 4)]
            colnames(data_bed[[i]]) = c("chr", "start", "end" , "name")
            data_bed[[i]]$strand = "+"
        } else{
            ##we have bed-3
            colnames(data_bed[[i]]) = c("chr", "start", "end")
            data_bed[[i]]$name = paste(data_bed[[i]]$chr,
                                       data_bed[[i]]$start,
                                       data_bed[[i]]$end,
                                       sep = "_")
            data_bed[[i]]$strand = "+"
        }
        index = which(x == names(data_bed[i]))
        data_bed[[i]]$site = names(c)[index]
    }
    return(data_bed)
}

load_vcf_files <- function(vcf_files, c, x) {
    data_vcf <- lapply(vcf_files, read.table, stringsAsFactors = FALSE)
    names(data_vcf) <- vcf_files

    for (i in seq_along(data_vcf)) {
        data_vcf[[i]] = data_vcf[[i]][, c(1, 2, 3, 4, 5)]
        colnames(data_vcf[[i]]) = c("chr", "end", "name", "ref", "alt")
        data_vcf[[i]]$start = (data_vcf[[i]]$end) - 1
        data_vcf[[i]]$strand = "+"
        index = which(x == names(data_vcf[i]))
        data_vcf[[i]]$site = names(c)[index]
        ##calculate the real end from the ALT col of vcf files
        data_vcf[[i]] = detect_vcf_end(data_vcf[[i]])
        ##reorder cols
        data_vcf[[i]] = data_vcf[[i]][, c("chr", "start", "end",
                                          "name", "strand", "site")]
    }
    return(data_vcf)
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

detect_vcf_end <- function(x) {
    indices = which((grepl(",", x$alt, fixed = TRUE)) == TRUE)
    index = which(nchar(x$alt) > 1)
    index = index[which(!(index %in% indices))]
    if (length(indices) == 0) {
        return(x)
    }

    for (i in 1:length(indices)) {
        s = (x$alt[indices[i]])
        alt = unlist(strsplit(s, ","))
        max = max(nchar(alt))
        x$end[indices[i]] = x$end[indices[i]] + max - 1
    }
    for (i in 1:length(index)) {
        x$end[index[i]] = x$end[index[i]] + nchar(x$start[index[i]]) - 1
    }
    return(x)
}



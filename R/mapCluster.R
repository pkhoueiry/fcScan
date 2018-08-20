#' @title fcScan

#' @export mapCluster

#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

mapCluster <- function(x, kent = NULL, chrom.sizes , chain) {
    ##x : clusters
    ##chain : path to chain file
    ##kent : path of kentutils
    ##chrom.sizes : path of chrom.sizes file
    ##assuming we have a GRange object as the final result

    input.bed <- data.frame(chrom = x[,1], start = x[,2],
                            end = x[,3], name = x[,4], score = x[,5], strand = x[,6] )
    tmp_dir = tempdir()
    write.table(input.bed,paste0(tmp_dir,"/input.bed"),
                sep = "\t",quote = FALSE, row.names = FALSE,col.names  = FALSE)

    ##apply bedToPsl on our indexed input files in order to
    ## convert each of them to a .psl file
    ##apply PSLMAP on each input.bed

    if(is.null(kent)){
        system(sprintf("bedToPsl %s %s %s",chrom.sizes ,
                       paste0(tmp_dir,"/input.bed"),paste0(tmp_dir,"/input.psl")))

        system(sprintf("pslMap -chainMapFile -swapMap -mapInfo=%s %s %s %s",
                       paste0(tmp_dir,"/results.txt"), paste0(tmp_dir,"/input.psl"),
                       chain ,paste0(tmp_dir,"/output.psl")))
    } else{
        system(sprintf("%s/bedToPsl %s %s %s",kent,chrom.sizes,
                       paste0(tmp_dir,"/input.bed"),paste0(tmp_dir,"/input.psl")))
        system(sprintf("%s/pslMap -chainMapFile -swapMap -mapInfo=%s %s %s %s",
                       kent, paste0(tmp_dir,"/results.txt"), paste0(tmp_dir,"/input.psl"),chain
                       ,paste0(tmp_dir,"/output.psl")))
    }

    results = read.delim(paste0(tmp_dir,"/results.txt"),
                         stringsAsFactors = FALSE ,fill = TRUE)
    ##mappedTName mappedTStart mappedTEnd mappedStrand mappedAligned
    ##srcTName srcTStart srcTEnd srcStrand
    results = results[,c(21,22,23,24,5,6,7,8)]
    return(results)
}




cache.getOrDownload <- function(url, rname) {
    bfc <- BiocFileCache(ask = FALSE)

    # check if url is tracked in local cache
    result <- bfcquery(bfc, url, field = "fpath")

    if (bfccount(result) == 0L) {

        # download the file located at url and add to cache
        cat("Downloading Reference Distribution: ")
        cachedFilePath <- bfcadd(bfc, rname = rname, fpath = url)
        
    } else {

        # get the ID if the file is already in the cache
        rid <- result[["rid"]]
        cachedFilePath <- bfcpath(bfc, rid)

        # check if the cached file has expired or needs update
        # (procedure recommended in the BiocFileCache Vignette)
        update_required <- bfcneedsupdate(bfc, rid)
        # update_required can be NA if it cannot be determined
        if (is.na(update_required)) {
            # TODO: while hosting on github, timestamps cannot be determined
            # correctly. Ignore uncorrectly determined update status for now
            update_required <- FALSE
        }
        if (update_required){
            cat("Updating Reference Distribution: ")
            cachedFilePath <- bfcdownload(bfc, rid, ask = FALSE)
        }
    }
    
    return(cachedFilePath)
    
}

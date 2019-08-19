
.cache.getOrDownload <- function(url, rname) {
    bfc <- BiocFileCache(ask = FALSE)

    # check if url is tracked in local cache
    result <- bfcquery(bfc, url, field="fpath", exact=TRUE)

    if (bfccount(result) == 0L) {

        # download the file located at url and add to cache
        print("Downloading Reference Distribution: \n")
        cachedFilePath <- bfcadd(bfc, rname = rname, fpath = url)

    } else {

        # get the ID if the file is already in the cache
        rid <- result[["rid"]]
        cachedFilePath <- bfcpath(bfc, rid)


        # TODO: Find hosting where the bfcneedsupdate functionality can
        # be utilized. While hosting this file on github, timestamps 
        # and expiry dates cannot be determined correctly. 
        # We'll ignore for now ..
        update_required <- FALSE

        # check if the cached file has expired or needs update
        # (procedure recommended in the BiocFileCache Vignette)
        # update_required <- bfcneedsupdate(bfc, rid)
        # update_required can be NA if it cannot be determined
        if (is.na(update_required)) {
            update_required <- FALSE
        }
        if (update_required){
            print("Updating Reference Distribution: \n")
            cachedFilePath <- bfcdownload(bfc, rid, ask = FALSE)
        }
    }
    return(cachedFilePath)
}

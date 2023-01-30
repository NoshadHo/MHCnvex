.HLACNV.addLogRatio <- function(tile, HLACNV.opts, opts) {
    if(HLACNV.opts$lr.mode == "pair" & any(!is.na(tile$n.cov))) {
        # tile$lr <- .HLACNV.calcLR(tile$t.cov, tile$n.cov)
        tile$lr <- cnvex:::.polishLogRatio(.HLACNV.calcLR(tile$t.cov, tile$n.cov), tile, opts)
    } else if(HLACNV.opts$lr.mode == "mean") {
        # tile$lr <- .HLACNV.calcLR(tile$t.cov, 1)
        tile$lr <-cnvex:::.polishLogRatio(.HLACNV.calcLR(tile$t.cov, 1), tile, opts)
    }

    if(!("lr" %in% colnames(mcols(tile)))) {
        stop("Lr mode is not correct, choose from 'mean' or 'pair'\n
             only select 'pair' if normal sample available")
    }
    return(tile)

}


.HLACNV.calcLR <- function(t.cov, n.cov) {
    lr.raw <- log2(t.cov / n.cov + 0.001)
    lr.raw <- ifelse(is.finite(lr.raw), lr.raw, NA_real_)
    return(lr.raw)
}

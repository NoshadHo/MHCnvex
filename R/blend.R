.HLACNV.blend.cnvTile.hla <- function(exons, cnv, HLACNV.opts) {
    ## exons columns needs for blending
    exons.toblend <- data.table(seqnames = exons$hla.type, start = 1:nrow(exons), end = (1:nrow(exons))+0.5, width=".", strand = "*",
                                n.cov.raw = exons$hap1.n.cov.raw+exons$hap2.n.cov.raw,
                                t.cov.raw = exons$hap1.t.cov.raw+exons$hap2.t.cov.raw,
                                n.cov = exons$hap1.n.cov+exons$hap2.n.cov,
                                t.cov = exons$hap1.t.cov+exons$hap2.t.cov,
                                gc = exons$gc,
                                target = TRUE, hq = TRUE, group = "HLA", arm = "chr6p", unmasked = 1, gap=0, blacklist=0, nC = 2, peak.id = exons$peak.id)

    ## cnv$tile columns needs for blending
    cnv.toblend <- copy(cnv$tile[,c("n.cov.raw","t.cov.raw","n.cov", "t.cov", "gc","target", "hq",
                                    "arm", "unmasked", "gap", "blacklist", "nC")])
    cnv.toblend$group   <- "GRCh38"
    cnv.toblend$peak.id <- NA
    cnv.toblend$idx     <- 1:length(cnv.toblend)

    ## overlaping the cnv$tile with hla exons
    coords <- .HLACNV.getConsts()
    cnv.coords.hits <- findOverlaps(cnv.toblend, coords, select = "first")

    ## finding exons idx to be blended based on GRCh38 coords
    coords.exons <- lapply(1:length(coords), function(i){
        row  <- coords[i]
        gene <- strsplit(row$gene_name, "-")[[1]][1]
        row$exon.tiles <- list(which(exons.toblend$seqnames == gene))
        return(row)
    })
    coords.exons <- do.call('c', coords.exons)

    ## let the blend begin ...
    blend <- .HLACNV.blend.inner(cnv.toblend, exons.toblend, cnv.coords.hits, coords.exons)

    return(blend)

}

.HLACNV.blend.cnvVar.hla <- function(mismatch, cnv, HLACNV.opts) {
    var.toblend <- data.table(seqnames = mismatch$hla.type, start = 1:nrow(mismatch), end = (1:nrow(mismatch))+0.5, width=".", strand = "*",
                              n.GT = "0/1", n.AF = mismatch$BAF.n, n.DP = mismatch$hap1.n.cov.raw+mismatch$hap2.n.cov.raw,
                              t.GT = "0/1", t.AF = mismatch$BAF.t, t.DP = mismatch$hap1.t.cov.raw+mismatch$hap2.t.cov.raw,
                              n.PASS = mismatch$baf.PASS, t.PASS = mismatch$hetero, SOURCE = cnv$var$SOURCE[1], TYPE = "SNP",
                              QUAL = max(cnv$var$QUAL), mask.loose=0, mask.strict=0, group = "HLA", peak.id = mismatch$peak.id)

    coords <- .HLACNV.getConsts()
    cnv.coords.hits <- findOverlaps(cnv$var, coords, select = "first")
    # cnv$tile

    coords.var <- lapply(1:length(coords), function(i){
        row  <- coords[i]
        gene <- strsplit(row$gene_name, "-")[[1]][1]
        row$exon.tiles <- list(which(var.toblend$seqnames == gene))
        return(row)
    })
    coords.var <- do.call('c', coords.var)

    if ("n.AF" %in% colnames(mcols(cnv$var))) { ## if normal available
        cnv.toblend <- (cnv$var[,c("n.GT","n.AF", "n.DP" ,
                                   "t.GT", "t.AF", "t.DP",
                                   "n.PASS", "t.PASS", "SOURCE", "TYPE", "QUAL", "mask.loose", "mask.strict")])
    } else {
        cnv.toblend <- (cnv$var[,c("t.GT", "t.AF", "t.DP",
                                   "t.PASS", "SOURCE", "TYPE", "QUAL", "mask.loose", "mask.strict")])
    }

    cnv.toblend$group    <- "GRCh38"
    cnv.toblend$peak.id  <- NA
    cnv.toblend$idx      <- 1:length(cnv.toblend)

    ## let the blend begin ...
    blend <- .HLACNV.blend.inner(cnv.toblend, var.toblend, cnv.coords.hits, coords.var)
    if (any(!is.na(blend$n.AF))) { ## if normal available
        blend$PASS <- blend$n.PASS & blend$t.PASS
    } else {
        blend$PASS <- blend$t.PASS
    }

    return(blend)
}

.HLACNV.blend.inner <- function(cnv.toblend, hla.toblend, cnv.coords.hits, coords.hla) {
    for(hla.idx in 1:length(coords.hla)) {
        row         <- coords.hla[hla.idx]
        cnv.idx     <- which(cnv.coords.hits == hla.idx)
        exon.idx    <- row$exon.tiles[[1]]
        exon.tiles  <- makeGRangesFromDataFrame(hla.toblend[exon.idx], keep.extra.columns = TRUE)

        if (length(cnv.idx) == 0) { ## TODO: make sure this doesn't cause any problems
            hits.idxs     <- unique(cnv.coords.hits[!is.na(cnv.coords.hits)])
            hla.idx.pos.n <- suppressWarnings(hits.idxs[min(which(hla.idx < hits.idxs))])   # what is next id
            hla.idx.pos.p <- suppressWarnings(hits.idxs[max(which(hla.idx > hits.idxs))])   # what is previous id
            ## find the tile in between:
            cnv.idx       <- suppressWarnings(fcase(!is.na(hla.idx.pos.n), as.integer(min(which(cnv.coords.hits == hla.idx.pos.n))+(hla.idx - hla.idx.pos.n)),
                                                    is.na(hla.idx.pos.n)&(!is.na(hla.idx.pos.p)), as.integer(max(which(cnv.coords.hits == hla.idx.pos.p))+(hla.idx - hla.idx.pos.p))))
        }

        if (length(exon.tiles) > 0) { ## TODO: make sure this doesn't cause any problems
            ## if the length is 0, then the idx is not important cause there is nothing to be added
            exon.tiles$idx <- min(cnv.idx)+0.0001*(1:length(exon.tiles))
        }

        cnv.toblend[cnv.toblend$idx %in% cnv.idx] <- NULL
        cnv.toblend <- suppressWarnings(c(cnv.toblend[cnv.toblend$idx < (min(cnv.idx)+1)],
                                          exon.tiles, cnv.toblend[cnv.toblend$idx > (max(cnv.idx)-1)]))
    }
    return(cnv.toblend)
}

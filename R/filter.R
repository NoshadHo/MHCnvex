.HLACNV.filterExons <- function(exons, HLACNV.opts) {
    exons <- copy(exons)

    filter.max.n.dev               <- HLACNV.opts$filter.perBase.max.n.dev
    filter.min.n.cov.sd            <- HLACNV.opts$filter.perBase.min.n.cov.sd
    filter.min.read.num            <- HLACNV.opts$filter.perBase.min.read.num
    filter.min.gene.hetero.percent <- HLACNV.opts$filter.perBase.min.gene.hetero.percent
    filter.min.t.cov.sd            <- HLACNV.opts$filter.perBase.min.t.cov.sd
    # normal baf
    if(any(!is.na(exons$BAF.n))) { ## if normal sample available
        exons[, baf.PASS := abs(BAF.n - 0.5) < filter.max.n.dev]
        exons[, cov.PASS := (hap1.n.cov.raw+hap2.n.cov.raw > mean(hap1.n.cov.raw+hap2.n.cov.raw, na.rm=T)-
                                 sd(hap1.n.cov.raw+hap2.n.cov.raw, na.rm=T)*filter.min.n.cov.sd),
              by= hla.type]
        exons[(hap1.n.cov.raw+hap2.n.cov.raw) < filter.min.read.num, cov.PASS := FALSE]
        exons[(hap1.t.cov.raw+hap2.t.cov.raw) < filter.min.read.num, cov.PASS := FALSE]
        exons[, hetero := ifelse(mean((baf.PASS[cov.PASS==TRUE]),na.rm=T) < filter.min.gene.hetero.percent,FALSE,TRUE), by=.(hla.type)]
    } else {
        exons[, baf.PASS := TRUE]
        exons[, cov.PASS := (hap1.t.cov+hap2.t.cov > mean(hap1.t.cov+hap2.t.cov, na.rm=T)-
                                 sd(hap1.t.cov+hap2.t.cov, na.rm=T)*filter.min.t.cov.sd),
              by= hla.type]
        exons[, hetero := TRUE]
    }
    exons[is.na(hetero), hetero := FALSE]

    return(exons)
}

.HLACNV.filterPeaks <- function(peaks, cov.tbl, HLACNV.opts) {
    ## filtering
    cov.tbl.split <- split(cov.tbl, cov.tbl$hla.type)
    peaks.split   <- split(peaks, peaks$gene)
    peaks.cov <- lapply(1:length(peaks.split), function(i) {

        subpeaks     <- peaks.split[[i]][peak.cov.PASS == TRUE]
        gene.covs    <- cov.tbl.split[[subpeaks$gene[1]]]
        ## finding adjacent regions of each peak
        start.idx <- subpeaks$peak.index - HLACNV.opts$filter.perBase.peak.adjacent
        start.idx <- ifelse(start.idx < 1, 1, start.idx) # if we the requested range not available

        end.idx   <- subpeaks$peak.index + HLACNV.opts$filter.perBase.peak.adjacent
        end.idx   <- ifelse(end.idx > nrow(gene.covs), nrow(gene.covs), end.idx) # if we the requested range not available
        ranges    <- data.table(start = gene.covs[start.idx]$start, end = gene.covs[end.idx]$start)
        setkey(ranges, start, end)

        overlap <- foverlaps(gene.covs[, end := start], ranges, by.x = c("start", "end"), by.y = c("start", "end"),which = TRUE, mult = "all",nomatch = NULL)

        gene.covs[overlap$xid, peak.adj.PASS := TRUE]
        gene.covs[is.na(peak.adj.PASS), peak.adj.PASS := FALSE]
        gene.covs[overlap$xid, peak.id := subpeaks$id[overlap$yid]]

        return(gene.covs)
    })
    peaks.cov <- rbindlist(peaks.cov)

    return(peaks.cov)
}


.HLACNV.MHCblacklist <- function(cnv, opts, HLACNV.opts) {
    ## cnvex ordinary blacklists without giab
    cnv$tile$hq <- TRUE
    t.cov.raw <- cnv$tile$t.cov.raw
    t.cov.raw[is.na(t.cov.raw)] <- 0
    n.cov.raw <- cnv$tile$n.cov.raw
    n.cov.raw[is.na(n.cov.raw)] <- 0
    tile.totalcov <- (t.cov.raw + n.cov.raw) * width(cnv$tile)
    cnv$tile$hq <-
        cnv$tile$gap <  opts$tile.hq.max.gap &
        cnv$tile$unmasked >  opts$tile.hq.min.unmasked &
        cnv$tile$blacklist <  opts$tile.hq.max.blacklist &
        # cnv$tile$giab.difficults < opts$tile.hq.max.giab.difficults & ## TODO: this threshold needs tuning
        tile.totalcov > opts$tile.hq.min.totalcov

    # cnv$tile$target <- ifelse(cnv$tile$arm %in% c("chr6p", "chr6q"), cnv$tile$target, FALSE)

    ## add MHCnvex specific blacklist
    if (HLACNV.opts$filter.MHC.blacklist == TRUE) {
        MHCblacklist <- system.file("extdata/GRCh38/mhc_blacklist.bed", package="MHCnvex")
        MHCblacklist <- cnvex:::.robust.import(MHCblacklist, seqinfo(cnv$tile))
        MHCblacklist <- GenomicRanges::reduce((MHCblacklist))

        tmp <- findOverlaps(cnv$tile, MHCblacklist)
        tmp <- data.table(
            tile=queryHits(tmp),
            blacklist=width(pintersect(cnv$tile[queryHits(tmp)], MHCblacklist[subjectHits(tmp)]))
        )
        setkey(tmp, tile)
        tmp <- tmp[J(seq_along(cnv$tile))]
        tmp[is.na(blacklist), blacklist:=0]
        tmp <- tmp[,.(blacklist=sum(blacklist)),by=tile]
        cnv$tile$MHCblacklist <- tmp$blacklist / width(cnv$tile)

        ## combine
        cnv$tile$hq <- cnv$tile$hq & (cnv$tile$MHCblacklist < opts$tile.hq.max.blacklist)
    }

    return(cnv)
}

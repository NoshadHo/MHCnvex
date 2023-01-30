.HLACNV.makeExonsPerBase <- function(basepair, ref, base.lr.mode, HLACNV.opts, peak.alignid=NULL) {

    align.basepair <- .HLACNV.basepairMigrate2Align(basepair, ref)

    exons      <- .HLACNV.basepair2exon(align.basepair, HLACNV.opts)
    exons      <- .HLACNV.addAF(exons, HLACNV.opts)

    exons      <- .HLACNV.filterExons(exons, HLACNV.opts)
    # exons      <- .HLACNV.addPeak(exons,HLACNV.opts)

    ## if pd provided, use pool peaks
    exons      <- .HLACNV.addPeak(exons,HLACNV.opts, peak.alignid)
    
    if(base.lr.mode == "peakTile") {
        ## adding gc per peak
        peak.gc <- .HLACNV.calcGcPerBaseMethod(exons)

        exons <- .HLACNV.basepair2Peak(exons, HLACNV.opts, align.basepair)

        exons <- merge.data.table(exons, peak.gc, by = "peak.id", all.x = TRUE)
    }

    if(HLACNV.opts$exons.fix.homo == TRUE) {
        exons <- .HLACNV.fixForHomozygous(exons, HLACNV.opts)
    }

    ## make sure missed peaks not causing any problem
    if (!all(is.na(exons$peak.id))) {
        exons[tstrsplit(peak.id,"-")[[2]] == "7410CCC", ':='(baf.PASS = FALSE, cov.PASS = FALSE, peak.PASS = FALSE)]
    } else {
        warning("No Peak were found!!!\nChange the minimum coverage setting")
    }

    ##
    exons$align.id <- exons$start
    return(exons)

}

.HLACNV.basepairMigrate2Align <- function(basepair, ref) {
    ## base.lr.mode: peakTile or peakBase
    align.basepair  <- copy(basepair)

    ## first migrate to align coords
    mapping.tbl <- ref$mapping.tbl


    align.basepair <- merge.data.table(align.basepair, mapping.tbl[, c("seqnames", "align.id")], all = FALSE,
                                   by.x = c("seqnames", "align.id"), by.y = c("seqnames", "align.id")) #start=end
    align.basepair <- align.basepair[!is.na(align.id)] ## remove the ones that are in raw but not in imp
    align.basepair$start      <- align.basepair$align.id
    align.basepair$end        <- align.basepair$align.id

    return(align.basepair[,.SD, .SDcols = colnames(basepair)])
}

.HLACNV.basepair2exon <- function(align.basepair, HLACNV.opts) {
    basepair.type <- split(align.basepair, align.basepair$hla.type)
    exon <- lapply(basepair.type, function(subbasepair) {
        alleles     <- split(subbasepair, subbasepair$seqnames)
        allele1     <- alleles[[1]]
        allele2     <- alleles[[2]]
        common.pos  <- base::intersect(allele1$align.id, allele2$align.id)
        allele1     <- allele1[allele1$align.id %in% common.pos]
        allele2     <- allele2[allele2$align.id %in% common.pos]

        return(data.table(hap1=allele1$seqnames, hap2=allele2$seqnames, seq1 = allele1$seq, seq2 = allele2$seq, hla.type=allele1$hla.type,
                          align.id = allele1$align.id,
                          start = allele1$start, end = allele2$end,
                          hap1.n.cov.raw = allele1$n.cov.raw, hap2.n.cov.raw = allele2$n.cov.raw,
                          hap1.t.cov.raw = allele1$t.cov.raw, hap2.t.cov.raw = allele2$t.cov.raw,
                          hap1.n.cov = allele1$n.cov, hap2.n.cov = allele2$n.cov,
                          hap1.t.cov = allele1$t.cov, hap2.t.cov = allele2$t.cov,
                          is.exon=allele1$is.exon, ref.exon.id = allele1$ref.exon.id))
    })
    exon <- rbindlist(exon)

    return(exon)
}


.HLACNV.addAF <- function(exons, HLACNV.opts) {
    exons$BAF.n <- exons$hap1.n.cov/ (exons$hap1.n.cov + exons$hap2.n.cov + 0.001) # 0.001 to prevent 0/0
    exons$BAF.t <- exons$hap1.t.cov/ (exons$hap1.t.cov + exons$hap2.t.cov + 0.001) # 0.001 to prevent 0/0

    return(exons)
}

.HLACNV.makeBAFPerBase <- function(align.basepair, HLACNV.opts) {
    basepair.type <- split(align.basepair, align.basepair$hla.type)
    BAF <- lapply(basepair.type, function(subbasepair) {
        alleles     <- split(subbasepair, subbasepair$seqnames)
        allele1     <- alleles[[1]]
        allele2     <- alleles[[2]]
        common.pos  <- base::intersect(allele1$align.id, allele2$align.id)
        allele1     <- allele1[allele1$align.id %in% common.pos]
        allele2     <- allele2[allele2$align.id %in% common.pos]

        BAF.n <- allele1$n.cov/ (allele1$n.cov + allele2$n.cov + 0.001) # 0.001 to prevent 0/0
        BAF.t <- allele1$t.cov/ (allele1$t.cov + allele2$t.cov + 0.001) # 0.001 to prevent 0/0
        gc    <- NA
        return(data.table(hap1=allele1$seqnames, hap2=allele2$seqnames, seq1 = allele1$seq, seq2 = allele2$seq, hla.type=allele1$hla.type,
                          start = allele1$start, end = allele2$end,
                          hap1.n.cov.raw = allele1$n.cov.raw, hap2.n.cov.raw = allele2$n.cov.raw,
                          hap1.t.cov.raw = allele1$t.cov.raw, hap2.t.cov.raw = allele2$t.cov.raw,
                          hap1.n.cov = allele1$n.cov, hap2.n.cov = allele2$n.cov,
                          hap1.t.cov = allele1$t.cov, hap2.t.cov = allele2$t.cov,
                          n.baf=BAF.n, t.baf=BAF.t, gc = gc,
                          is.exon=allele1$is.exon, ref.exon.id = allele1$ref.exon.id))
    })
    BAF <- rbindlist(BAF)
    return(BAF)

}

.HLACNV.makeLrPerBase <- function(align.basepair, HLACNV.opts) {
    basepair.type <- split(align.basepair, align.basepair$hla.type)
    LR <- lapply(basepair.type, function(subbasepair) {
        alleles     <- split(subbasepair, subbasepair$seqnames)
        allele1     <- alleles[[1]]
        allele2     <- alleles[[2]]
        common.pos  <- base::intersect(allele1$align.id, allele2$align.id)
        allele1     <- allele1[allele1$align.id %in% common.pos]
        allele2     <- allele2[allele2$align.id %in% common.pos]

        lr.dip <- log2((allele1$t.cov + allele2$t.cov)/ (allele1$n.cov + allele2$n.cov + 0.001))
        lr.hap1 <- log2((allele1$t.cov)/ (allele1$n.cov + 0.001))
        lr.hap2 <- log2((allele2$t.cov)/ (allele2$n.cov + 0.001))
        return(data.table(hap1=allele1$seqnames, hap2=allele2$seqnames, seq1 = allele1$seq, seq2 = allele2$seq, hla.type=allele1$hla.type, start = allele1$start, end = allele2$end,
                          hap1.n.cov = allele1$n.cov, hap2.n.cov = allele2$n.cov,
                          hap1.t.cov = allele1$t.cov, hap2.t.cov = allele2$t.cov,
                          lr.dip = lr.dip, lr.hap1 = lr.hap1, lr.hap2 = lr.hap2))
    })
    LR <- rbindlist(LR)

    return(LR)
}

.HLACNV.fixForHomozygous <- function(exons, HLACNV.opts) {
    ## fix the cov
    exons[hetero==FALSE, hap1.n.cov := hap1.n.cov+hap2.n.cov]
    exons[hetero==FALSE, hap1.t.cov := hap1.t.cov+hap2.t.cov]

    exons[hetero==FALSE, hap2.n.cov := 0]
    exons[hetero==FALSE, hap2.t.cov := 0]

    ## fix the baf filter
    exons[hetero==FALSE, baf.PASS := FALSE]

    ## fix the baf
    if ("n.baf" %in% colnames(exons)) {
        exons[hetero==FALSE, n.baf := 1]
    }
    if ("t.baf" %in% colnames(exons)) {
        exons[hetero==FALSE, t.baf := 1]
    }

    ## fix for LR
    if ("lr.dip" %in% colnames(exons)) {
        exons[hetero==FALSE, lr.dip := log2((hap1.t.cov + hap2.t.cov)/ (hap1.n.cov + hap2.n.cov))]
    }
    if ("lr.hap1" %in% colnames(exons)) {
        exons[hetero==FALSE, lr.hap1 := log2((hap1.t.cov)/ (hap1.n.cov + hap2.n.cov))]
    }
    if ("lr.hap2" %in% colnames(exons)) {
        exons[hetero==FALSE, lr.hap2 := log2((hap2.t.cov)/ (hap1.n.cov + hap2.n.cov))]
    }

    return(exons)

    return(exons)
}

.HLACNV.makeBAF <- function(basepair, HLACNV.opts) {
    basepair.type <- split(basepair, basepair$hla.type)
    BAF <- lapply(basepair.type, function(subbasepair) {
        alleles <- split(subbasepair, subbasepair$seqnames)
        BAF.n <- alleles[[1]]$n.cov/ (alleles[[1]]$n.cov + alleles[[2]]$n.cov)
        BAF.t <- alleles[[1]]$t.cov/ (alleles[[1]]$t.cov + alleles[[2]]$t.cov)
        gc    <- (alleles[[1]]$gc + alleles[[2]]$gc)/2
        return(data.table(hap1=alleles[[1]]$seqnames, hap2=alleles[[2]]$seqnames, hla.type=alleles[[1]]$hla.type,
                          hap1.n.cov = alleles[[1]]$n.cov, hap2.n.cov = alleles[[2]]$n.cov,
                          hap1.t.cov = alleles[[1]]$t.cov, hap2.t.cov = alleles[[2]]$t.cov,
                          n.baf=BAF.n, t.baf=BAF.t, exon.idx=alleles[[1]]$idx, gc = gc))
    })
    BAF <- rbindlist(BAF)
    return(BAF)
}

.HLACNV.makeLr <- function(basepair,HLACNV.opts) {
    basepair.type <- split(basepair, basepair$hla.type)
    LR <- lapply(basepair.type, function(subbasepair) {
        alleles <- split(subtile, subbasepair$seqnames)
        lr.dip <- log2((alleles[[1]]$t.cov + alleles[[2]]$t.cov)/ (alleles[[1]]$n.cov + alleles[[2]]$n.cov))
        lr.hap1 <- log2((alleles[[1]]$t.cov)/ (alleles[[1]]$n.cov))
        lr.hap2 <- log2((alleles[[2]]$t.cov)/ (alleles[[2]]$n.cov))
        return(data.table(hap1=alleles[[1]]$seqnames, hap2=alleles[[2]]$seqnames, hla.type=alleles[[1]]$hla.type,
                          hap1.n.cov = alleles[[1]]$n.cov, hap2.n.cov = alleles[[2]]$n.cov,
                          hap1.t.cov = alleles[[1]]$t.cov, hap2.t.cov = alleles[[2]]$t.cov,
                          lr.dip = lr.dip, lr.hap1 = lr.hap1, lr.hap2 = lr.hap2,
                          exon.idx=alleles[[1]]$idx))
    })
    LR <- rbindlist(LR)
}

.HLACNV.calcLR <- function(t.cov, n.cov) {
    lr.raw <- log2(t.cov / n.cov+ 0.001)
    lr.raw <- ifelse(is.finite(lr.raw), lr.raw, NA_real_)
    return(lr.raw)
}

.HLACNV.calcGcPerBaseMethod <- function(exons) {
    exon.peaks <- split(exons, exons$peak.id)
    gc.content <- lapply(exon.peaks, function(subexon.peaks) {
        seq1    <- Biostrings::DNAString(paste(subexon.peaks$seq1, collapse = ""))
        seq1.gc <- letterFrequency(seq1, "GC", as.prob = TRUE)
        seq2    <- Biostrings::DNAString(paste(subexon.peaks$seq2, collapse = ""))
        seq2.gc <- letterFrequency(seq2, "GC", as.prob = TRUE)

        seq.gc <- mean(seq1.gc, seq2.gc)

        return(seq.gc)
    })
    return(data.table(peak.id = names(gc.content), gc = unlist(gc.content)))
}

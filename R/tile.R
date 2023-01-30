.HLACNV.makeBasePair <- function(ref) {
    basepair <- lapply(1:nrow(ref$seq.imputed), function(i) {
        row <- ref$mapping.tbl[seqnames == ref$seq.imputed[i]$name]

        dt  <- data.table(seqnames = row$seqnames[1], start = row$imp.id, end = row$imp.id, seq=row$seq.raw, align.id = row$align.id,
                          hla.type = row$hla.type, is.exon = row$is.exon, exon.start = row$exon.start, exon.end = row$exon.start,
                          ref.exon.id = row$ref.exon.id)
        return(dt)
    })
    basepair <- rbindlist(basepair)
    # basepair$hla.type <- do.call(rbind,strsplit(basepair$seqnames, "\\*"))[,1]
    return(basepair)

}


.HLACNV.importBam <- function(t.bam, n.bam, basepair, fasta, HLACNV.opts) {
    ## normalize by sequencing depth
    if (!is.null(t.bam)) {
        t.cov.raw <- .HLACNV.runMosdepthBasepair(t.bam, basepair, fasta, HLACNV.opts$cores, HLACNV.opts$basepair.mapq, HLACNV.opts$basepair.cov.flags)
    } else {
        t.cov.raw <- NA_real_
    }
    if (!is.null(n.bam)) {
        n.cov.raw <- .HLACNV.runMosdepthBasepair(n.bam, basepair, fasta, HLACNV.opts$cores, HLACNV.opts$basepair.mapq, HLACNV.opts$basepair.cov.flags)
    } else {
        n.cov.raw <- NA_real_
    }

    basepair[, ':=' (
        n.cov.raw = n.cov.raw,
        t.cov.raw = t.cov.raw

    )]

    # basepair$t.cov.raw[is.na(basepair$t.cov.raw)] <- 0
    # basepair$n.cov.raw[is.na(basepair$n.cov.raw)] <- 0

    return(basepair)
}

.HLACNV.runMosdepthBasepair <- function(bam, basepair, fasta, cores, mapq, basepair.cov.flags) {
    bed <- tempfile("mosdepth_", fileext=".bed")
    rtracklayer::export(basepair, bed)
    out <- .HLACNV.runMosdepth(bam, bed, fasta, cores, mapq, basepair.cov.flags)
    unlink(bed)
    return(out)
}

.HLACNV.runMosdepth <- function(bam, by, fasta, cores, mapq, basepair.cov.flags) {
    if (is.character(by) && file.exists(by)) {
        cov.col <- "V5"
    } else if (!is.na(suppressWarnings(as.integer(by))) ) {
        cov.col <- "V4"
    } else {
        stop("by should be a file or window-size")
    }

    prefix <- tempfile("mosdepth_")
    ret <- system2("mosdepth", sprintf("-f %s -b %s -F %s -Q %s -n -t%s  %s %s", fasta, by, basepair.cov.flags, mapq, cores, prefix, bam))
    out.fn <- list.files(dirname(prefix), paste0(basename(prefix), ".regions.bed.gz$"),
                         full.names = TRUE)
    if (!ret & file.exists(out.fn)) {
        out.data <- with(fread(out.fn), GRanges(V1, IRanges(V2, V3), cov=get(cov.col)))
        cov <- sort(sortSeqlevels(out.data))$cov
    } else {
        stop(sprintf("mosdepth run failed: %s", ret))
    }
    out.fns <- list.files(dirname(prefix), paste0(basename(prefix)), full.names = TRUE)
    rets <- sapply(out.fns, unlink)
    if (any(rets)) {
        stop("could not remove temp. files")
    }
    return(cov)
}


.HLACNV.normCoverageBasepair <- function(basepair, t.normalization.param, n.normalization.param) {
    n.cov <- .HLACNV.normCoverage(basepair$n.cov.raw, basepair, n.normalization.param)
    t.cov <- .HLACNV.normCoverage(basepair$t.cov.raw, basepair, t.normalization.param)

    basepair$n.cov <- n.cov
    basepair$t.cov <- t.cov

    return(basepair)
}

.HLACNV.normCoverage <- function(cov, basepair, normalization.param) {
    ## normalized coverage
    cov.n <- cov*(basepair$end - basepair$start + 1)

    cov.n <- cov.n/(normalization.param/1e6)
    cov.n <- 1000*cov.n/(basepair$end - basepair$start + 1)

    return(cov.n)
}

.HLACNV.annotatebasepair <- function(basepair, seq.stringset, HLACNV.opts) {
    ## HLA type
    basepair$hla.type <- do.call(rbind,strsplit(as.character(basepair$seqnames),"\\*"))[,1]

    ## GC
    if(HLACNV.opts$mode == "exon") {
        tmp <- getSeq(seq.stringset, makeGRangesFromDataFrame(basepair))
        basepair$gc <- letterFrequency(tmp, "GC", as.prob=TRUE)[,1]
    } else if (HLACNV.opts$mode == "base") {
        basepair$gc <- NA
    }

    return(basepair)
}

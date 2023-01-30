.HLACNV.readPoolCoverage <- function(n.bam, fasta, calls, id, cov.tot.n, cov.tot.t, HLACNV.opts) {
    ## id: patient name, same as patient value to hapster

    ## make ref
    ref <- .HLACNV.readRef(fasta, calls, HLACNV.opts)

    ## make basepair
    basepair <- .HLACNV.makeBasePair(ref)
    basepair <- .HLACNV.importBam(NULL, n.bam, basepair, fasta, HLACNV.opts)
    basepair <- .HLACNV.normCoverageBasepair(basepair, cov.tot.t, cov.tot.n)

    exons  <- .HLACNV.makeExonsPerBase(basepair, ref, "peakBase", HLACNV.opts)
    exons[, n.total.cov := hap1.n.cov + hap2.n.cov]
    exons$sample.id <- id
    return(exons[, c("hap1", "hap2", "hla.type", "hap1.n.cov", "hap2.n.cov", "n.total.cov",
                     "align.id", "sample.id", "is.exon", "ref.exon.id", "seq1", "hetero", "baf.PASS", "cov.PASS")])

}


.HLACNV.pickSharedBasepairs <- function(cov.cohort) {
    ## making sure all the selected basepairs
    ## are shared among all samples

    #expected number of each align.id

    n.sample <- length(unique(cov.cohort$sample.id))

    basepair.n   <- cov.cohort[, .(n = .N), by=.(hla.type, align.id)]
    if(any(basepair.n$n > n.sample)) {
        stop("More than 1 file with same patient id. Making pool fail!")
    }
    basepair.n   <- basepair.n[n == n.sample]

    setkey(cov.cohort, hla.type, align.id)
    cov.shared <- cov.cohort[basepair.n[,-"n"]]

    return(cov.shared)
}


.HLACNV.ReadCohortCov <- function(hapster.fn, cnv.fn, HLACNV.opts, opts) {
    ## read cohort coverage
    cov.cohort <- foreach(i = 1:length(hapster.fn)) %dopar% {
        hapster <- hapster.fn[i]
        cnv <- cnv.fn[i]
        ## finding normal id
        sfiles <- list.files(hapster)
        n.id   <- sfiles[grepl("*.csv$", sfiles)]
        n.id   <- strsplit(n.id, "_haplotype.csv")[[1]]

        id    <- basename(hapster)
        bam.n <- sprintf("%s/alignments/%s_germline_imputed.bam",bam, n.id)
        fasta <- sprintf("%s/refs/%s_germline_imputed.fa",bam, id)
        calls <- sprintf("%s/calls",hapster)

        # cnv <- readRDS(sprintf("%s/%s.rds", cnv, id))
        cnv <- readRDS(cnv)
        cov.tot.n <- sum((cnv$tile$n.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])
        cov.tot.t <- sum((cnv$tile$t.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])

        cov.tbl <- .HLACNV.readPoolCoverage(bam.n, fasta, calls, id, cov.tot.n, cov.tot.t, HLACNV.opts)

        return(cov.tbl)
    }
    cov.cohort <- rbindlist(cov.cohort)
    cov.shared <- .HLACNV.pickSharedBasepairs(cov.cohort)
    cov.split <- split(cov.shared, cov.shared$sample.id)

    cov.gene <- cov.split[[1]]$hla.type
    cov.mx   <- do.call(rbind, lapply(cov.split, function(cov) {return(cov$n.total.cov)}))
    colnames(cov.mx) <- cov.split[[1]]$align.id

    return(list(cov.mx = cov.mx, cov.gene = cov.gene, align.id = cov.split[[1]]$align.id))

}

### allelewise functions
.HLACNV.readPoolCoverageAllele <- function(n.bam, fasta, calls, id, cov.tot.n, cov.tot.t, HLACNV.opts) {
    ## id: patient name, same as patient value to hapster

    ## make ref
    ref <- .HLACNV.readRef(fasta, calls, HLACNV.opts)

    ## make basepair
    basepair <- .HLACNV.makeBasePair(ref)
    basepair <- .HLACNV.importBam(NULL, n.bam, basepair, fasta, HLACNV.opts)
    basepair <- .HLACNV.normCoverageBasepair(basepair, cov.tot.t, cov.tot.n)
    basepair$sample_id <- id

    exons  <- .HLACNV.makeExonsPerBase(basepair, ref, "peakBase", HLACNV.opts)
    exons[, n.total.cov := hap1.n.cov + hap2.n.cov]

    exons.info <- exons[ ,.(hetero = hetero[1]) ,by = hla.type]
    exons.info$sample.id <- id
    return(list(basepair = basepair, exons.info = exons.info))

}

.HLACNV.pickSharedBasepairsAllele <- function(cov.cohort) {
    ## making sure all the selected basepairs
    ## are shared among all samples

    #expected number of each align.id

    n.allele <- length(unique(cov.cohort$sample_id))*2

    basepair.n   <- cov.cohort[, .(n = .N), by=.(hla.type, align.id)]
    if(any(basepair.n$n > n.allele)) {
        stop("More than 1 file with same patient id. Making pool fail!")

    }
    basepair.n   <- basepair.n[n == n.allele]

    setkey(cov.cohort, hla.type, align.id)
    cov.shared <- cov.cohort[basepair.n[,-"n"]]

    return(cov.shared)
}

.HLACNV.geneSpecificPool <- function(pool, gene) {
    gspool <- list()

    gspool$median.cov    <- pool$median.cov[pool$gene == gene]
    gspool$gene          <- pool$gene[pool$gene == gene]
    gspool$align.id      <- pool$align.id[pool$gene == gene]
    gspool$ref.exons     <- pool$ref.exons[pool$gene == gene]
    gspool$cohort.hetero <- pool$cohort.hetero[hla.type == gene]
    gspool$row.coverages <- pool$row.coverages
    gspool$peak.id       <- pool$peak.id[pool$gene == gene]
    gspool$coverages     <- pool$coverages[,pool$gene == gene]

    return(gspool)

}

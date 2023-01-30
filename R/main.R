#' @export
HLACNV.main <- function(sample_id, n.bam, t.bam, fasta, calls, cnv, purity, ploidy, output.fn, HLACNV.opts, poold=list()) {

    ## read data
    if (HLACNV.opts$panel.mode == "wxs") {
        opts <- cnvex::getOpts(system.file(sprintf("extdata/settings/%s.R", "exome-pair"), package="cnvex"))
    } else {
        opts <- cnvex::getOpts(system.file(sprintf("extdata/settings/%s.R", "genome-pair-15x"), package="cnvex"))
    }

    cnv <- readRDS(cnv)
    ##### filter
    cnv     <- .HLACNV.MHCblacklist(cnv,opts,HLACNV.opts)
    #####

    cov.tot.n <- sum((cnv$tile$n.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])
    cov.tot.t <- sum((cnv$tile$t.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])
    p <- purity
    P <- ploidy
    # HLACNV.opts <- .HLACNV.getOpts()
    ref <- .HLACNV.readRef(fasta, calls, HLACNV.opts)


    ## make tile
    basepair <- .HLACNV.makeBasePair(ref)
    basepair <- .HLACNV.importBam(t.bam, n.bam, basepair, ref$fasta, HLACNV.opts)
    basepair <- .HLACNV.normCoverageBasepair(basepair, cov.tot.t, cov.tot.n)


    ## make exons object
    # exons  <- .HLACNV.makeExonsPerBase(basepair, ref, "peakTile", HLACNV.opts)
    exons  <- .HLACNV.makeExonsPerBase(basepair, ref, "peakTile", HLACNV.opts, peak.alignid = poold$peaks)

    ## gc normalize
    if (HLACNV.opts$bias.gc.mode == "on") {
        opts$bias.gc.adjust.span.on = 0.5
        normGClr  <- .HLACNV.normGCLr(cnv, exons, HLACNV.opts, opts)
        exons    <- normGClr[[1]]
    }

    ## TODO: integrate HLACNV.opts with cnvex opts

    opts$lr.tumor <- HLACNV.opts$lr.mode
    ## make model HLA
    # mismatch  <- .HLACNV.makeMismatchPerBase(basepair, ref, HLACNV.opts)
    mismatch  <- .HLACNV.makeMismatchPerBase(basepair, ref, HLACNV.opts, poold$peaks)
    mhla      <- .HLACNV.modelHLA(exons, mismatch, basepair, cnv, HLACNV.opts, opts)
    mhla$tile <- .HLACNV.addLogRatio(mhla$tile,HLACNV.opts, opts)
    seg       <- .HLACNV.segmentation(mhla, opts) ## only chr6p segmentation
    fit       <- .HLACNV.modelFit(mhla, seg, purity, ploidy, opts)
    gene      <- .HLACNV.geneOut(fit,seg,mhla)
    haploids  <- .HLACNV.haploidCnv(basepair, gene, mismatch)


    p.adj = .HLACNV.plotLrBaf(mhla, gene, seg, purity, ploidy, plot.type = c("t.AF", "n.AF"), HLACNV.opts)+ggtitle(sample_id)
    # return(mhla, gene, haploids, p.adj)

    if (!file.exists(output.fn)) {
        system(sprintf("mkdir %s", output.fn))
    }
    ggsave(sprintf("%s/%s.png",output.fn, sample_id),plot = p.adj, device = "png", width = 500 ,height = 250, units = "mm")


    fwrite(gene, sprintf("%s/%s-gene.csv",output.fn, sample_id))
    fwrite(haploids, sprintf("%s/%s-haploids.csv",output.fn, sample_id))
    return(p.adj)
}

#' @export
HLACNV.makePdDiploid <- function(hapster.fn, cnv.fn, cores = 1, HLACNV.opts) {
  ## cnv.fn: path to cnv.rds file
  ## hapster.fn: path to hapster scratch
  doParallel::registerDoParallel(cores)
  ## read cohort coverage
  cov.cohort <- foreach(i = 1:length(hapster.fn)) %dopar% {
    hapster <- hapster.fn[i]
    cnv <- cnv.fn[i]
    ## finding normal id
    sfiles <- list.files(hapster)
    n.id   <- sfiles[grepl("*.csv$", sfiles)]
    n.id   <- strsplit(n.id, "_haplotype.csv")[[1]]

    id    <- basename(hapster)
    # bam.n <- sprintf("%s/alignments/%s_germline_imputed.bam",hapster, n.id)
    # fasta <- sprintf("%s/refs/%s_germline_imputed.fa",hapster, id)
    n.bam <- sprintf("%s/alignments/%s_germline_imputed.bam",hapster, n.id)
    fasta <- sprintf("%s/refs/%s_germline_imputed.fa",hapster, id)
    calls <- sprintf("%s/calls",hapster)

    # cnv <- readRDS(sprintf("%s/%s.rds", cnv, id))
    cnv <- readRDS(cnv)
    cov.tot.n <- sum((cnv$tile$n.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])
    cov.tot.t <- sum((cnv$tile$t.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])

    cov.tbl <- .HLACNV.readPoolCoverage(n.bam, fasta, calls, id, cov.tot.n, cov.tot.t, HLACNV.opts)

    return(cov.tbl)
  }
  cov.cohort <- rbindlist(cov.cohort)
  cov.shared <- .HLACNV.pickSharedBasepairs(cov.cohort)
  cov.split <- split(cov.shared, cov.shared$sample.id)

  cov.gene <- cov.split[[1]]$hla.type
  cov.mx   <- do.call(rbind, lapply(cov.split, function(cov) {return(cov$n.total.cov)}))
  colnames(cov.mx) <- cov.split[[1]]$align.id


  ### find the peak for each gene
  cov.mean <- colMeans(cov.mx, na.rm = TRUE)
  align.id <- cov.split[[1]]$align.id
  cov.mean.split <- split(cov.mean, cov.gene)
  align.id.split <- split(align.id, cov.gene)
  peaks <- lapply(1:length(cov.mean.split), function(j) {
    cov   <- cov.mean.split[[j]]
    id    <- align.id.split[[j]]
    gene  <- names(cov.mean.split)[j]
    peak  <- .HLACNV.detectPeaksInner(id, cov, HLACNV.opts = HLACNV.opts)
    peakx <- peak$x[which(peak$y.hat[peak$i] > HLACNV.opts$pool.min.cov.normalized)]
    loc   <- peak$i[which(peak$y.hat[peak$i] > HLACNV.opts$pool.min.cov.normalized)]
    return(data.table(gene = gene, peak.align.id = list(peakx), peak.loc = list(loc),cov.len = length(cov)))
  })
  peaks <- rbindlist(peaks)
  setkey(peaks, gene)
  peaks <- peaks[J(data.table(gene = names(cov.mean.split)))]
  peaks[is.na(cov.len), cov.len := lengths(cov.mean.split[gene])]
  ## defined the adjacent tiles of the peak
  selected.tiles <- lapply(1:nrow(peaks), function(i) {
    gene     <- peaks[i,]$gene
    peak.loc <- peaks[i,]$peak.loc[[1]]
    len      <- peaks[i,]$cov.len
    id       <- align.id.split[[gene]]
    tiles <- lapply(peak.loc, function(loc) {
      max.loc <- ifelse(loc+HLACNV.opts$filter.perBase.peak.adjacent < len, loc+HLACNV.opts$filter.perBase.peak.adjacent, len)
      min.loc <- ifelse(loc-HLACNV.opts$filter.perBase.peak.adjacent > 1, loc-HLACNV.opts$filter.perBase.peak.adjacent, 1)
      return(data.table(peak.id = paste0(gene,"-",loc), loc = seq(min.loc, max.loc, 1), align.id = id[min.loc:max.loc]))
    })
    tiles <- rbindlist(tiles)


    if(!("align.id" %in% colnames(tiles))) { #if no peak was found
      tiles <- data.table(align.id=id, peak.id = NA, loc = 1:length(id),gene = gene)
    } else {
      ## handle if peak regions overlap
      tiles <- tiles[, .SD[1],by=align.id]
      setkey(tiles, align.id)
      tiles <- tiles[J(data.table(align.id = id))]
      tiles$gene <- gene
    }

    return(tiles)

  })
  selected.tiles <- rbindlist(selected.tiles)

  ## summarize everything into a pool
  poold <- list(median.cov = cov.mean, gene = cov.split[[1]]$hla.type, align.id = cov.split[[1]]$align.id,
                ref.exons = cov.split[[1]]$ref.exon.id)
  poold$peak.id   <- selected.tiles$peak.id
  poold$coverages <- cov.mx
  poold$peaks     <- data.table(align.id = cov.split[[1]]$align.id, hla.type = cov.split[[1]]$hla.type, peak.id = selected.tiles$peak.id)
  # saveRDS(pool, "/mctp/users/noshadh/data/agilent.hla.pool.rds")
  doParallel::stopImplicitCluster()
  return(poold)

}


#' @export
HLACNV.makePdAllele <- function(hapster.fn, cnv.fn, cores = 1, HLACNV.opts) {
  ## cnv.fn: path to cnv.rds cov.tblfile
  ## hapster.fn: path to hapster scratch
  HLACNV.opts.allele <- HLACNV.opts
  HLACNV.opts.allele$filter.perBase.min.read.num <- HLACNV.opts.allele$filter.perBase.min.read.num/2
  doParallel::registerDoParallel(cores)
  ## read cohort coverage
  cohort.info <- foreach(i = 1:length(hapster.fn)) %dopar% {
    hapster <- hapster.fn[i]
    cnv <- cnv.fn[i]
    ## finding normal id
    sfiles <- list.files(hapster)
    n.id   <- sfiles[grepl("*.csv$", sfiles)]
    n.id   <- strsplit(n.id, "_haplotype.csv")[[1]]

    id    <- basename(hapster)
    # bam.n <- sprintf("%s/alignments/%s_germline_imputed.bam",hapster, n.id)
    # fasta <- sprintf("%s/refs/%s_germline_imputed.fa",hapster, id)
    n.bam <- sprintf("%s/alignments/%s_germline_imputed.bam",hapster, n.id)
    fasta <- sprintf("%s/refs/%s_germline_imputed.fa",hapster, id)
    calls <- sprintf("%s/calls",hapster)

    # cnv <- readRDS(sprintf("%s/%s.rds", cnv, id))
    cnv <- readRDS(cnv)
    cov.tot.n <- sum((cnv$tile$n.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])
    cov.tot.t <- sum((cnv$tile$t.cov.raw * width(cnv$tile))[cnv$tile$target == TRUE])

    sample_info <- .HLACNV.readPoolCoverageAllele(n.bam, fasta, calls, id, cov.tot.n, cov.tot.t, HLACNV.opts.allele)

    return(sample_info)
  }
  cov.cohort <- rbindlist(lapply(cohort.info, '[[',1))
  cov.shared <- .HLACNV.pickSharedBasepairsAllele(cov.cohort)

  cov.split <- split(cov.shared, paste0(cov.shared$sample_id, "-", cov.shared$seqnames))

  ## make sample table
  cohort.hetero <- rbindlist(lapply(cohort.info, '[[',2))

  cov.gene <- cov.split[[1]]$hla.type
  cov.mx   <- do.call(rbind, lapply(cov.split, function(cov) {return(cov$n.cov)}))
  colnames(cov.mx) <- cov.split[[1]]$align.id


  ### find the peak for each gene
  cov.mean <- colMeans(cov.mx, na.rm = TRUE)
  align.id <- cov.split[[1]]$align.id
  cov.mean.split <- split(cov.mean, cov.gene)
  align.id.split <- split(align.id, cov.gene)
  peaks <- lapply(1:length(cov.mean.split), function(j) {
    cov   <- cov.mean.split[[j]]
    id    <- align.id.split[[j]]
    gene  <- names(cov.mean.split)[j]
    peak  <- .HLACNV.detectPeaksInner(id, cov, HLACNV.opts = HLACNV.opts.allele)
    peakx <- peak$x[which(peak$y.hat[peak$i] > HLACNV.opts.allele$pool.min.cov.normalized/2)] ## move the /2 part into setting file if decided to go with allele-wise pool
    loc   <- peak$i[which(peak$y.hat[peak$i] > HLACNV.opts.allele$pool.min.cov.normalized/2)] ## move the /2 part into setting file if decided to go with allele-wise pool
    return(data.table(gene = gene, peak.align.id = list(peakx), peak.loc = list(loc),cov.len = length(cov)))
  })
  peaks <- rbindlist(peaks)
  setkey(peaks, gene)
  peaks <- peaks[J(data.table(gene = names(cov.mean.split)))]
  peaks[is.na(cov.len), cov.len := lengths(cov.mean.split[gene])]
  ## defined the adjacent tiles of the peak
  selected.tiles <- lapply(1:nrow(peaks), function(i) {
    gene     <- peaks[i,]$gene
    peak.loc <- peaks[i,]$peak.loc[[1]]
    len      <- peaks[i,]$cov.len
    id       <- align.id.split[[gene]]
    tiles <- lapply(peak.loc, function(loc) {
      max.loc <- ifelse(loc+HLACNV.opts.allele$filter.perBase.peak.adjacent < len, loc+HLACNV.opts.allele$filter.perBase.peak.adjacent, len)
      min.loc <- ifelse(loc-HLACNV.opts.allele$filter.perBase.peak.adjacent > 1, loc-HLACNV.opts.allele$filter.perBase.peak.adjacent, 1)
      return(data.table(peak.id = paste0(gene,"-",loc), loc = seq(min.loc, max.loc, 1), align.id = id[min.loc:max.loc]))
    })
    tiles <- rbindlist(tiles)


    if(!("align.id" %in% colnames(tiles))) { #if no peak was found
      tiles <- data.table(align.id=id, peak.id = NA, loc = 1:length(id),gene = gene)
    } else {
      ## handle if peak regions overlap
      tiles <- tiles[, .SD[1],by=align.id]
      setkey(tiles, align.id)
      tiles <- tiles[J(data.table(align.id = id))]
      tiles$gene <- gene
    }

    return(tiles)

  })
  selected.tiles <- rbindlist(selected.tiles)

  ## summarize everything into a pool
  poold <- list(median.cov = cov.mean, gene = cov.split[[1]]$hla.type, align.id = cov.split[[1]]$align.id,
                ref.exons = cov.split[[1]]$ref.exon.id, cohort.hetero = cohort.hetero, row.coverages = rownames(cov.mx))
  poold$peak.id   <- selected.tiles$peak.id
  poold$coverages <- cov.mx
  poold$peaks     <- data.table(align.id = poold$align.id, hla.type = poold$gene, peak.id = poold$peak.id)
  # saveRDS(pd, "/mctp/users/noshadh/data/agilent.hla.pd.allele.rds")
  doParallel::stopImplicitCluster()
  return(poold)

}



# pdf(file=sprintf("/mctp/users/noshadh/work/cnv-rna-prot/results/MHCnvex/mhcnvexpool.pdf"),width = 10,height = 5, paper = "a4")
# ((MHCnvex:::.HLACNV.plotPool(poold, HLACNV.opts)))
# dev.off()

#
#
# foreach(i = 1:nrow(brodie)) %dopar% {
#     digest = readRDS(brodie$digest[i])
#     cnv = readRDS(brodie$cnv[i])
#     #
#     sample_id = brodie$patient[i]
#     n.bam = brodie$bam.n[i]
#     t.bam = brodie$bam.t[i]
#     fasta = brodie$fasta[i]
#     calls = brodie$calls[i]
#     purity = digest$model.purity
#     ploidy = digest$model.ploidy
#     digest = brodie$digest[i]
#     cnv = brodie$cnv[i]
#
#     output.fn <- paste0("/mctp/share/users/noshadh/brodie/loh_samples/", sample_id)
#     system(sprintf("mkdir %s", output.fn))
#
#     HLACNV.main(sample_id, n.bam, t.bam, fasta, calls, cnv, purity, ploidy, output.fn, HLACNV.opts)
# }


#
# hapster.fn <- list.files("/mctp/users/noshadh/work/HLACNV/data/cohort_177/cohort_177_hapster", full.names = TRUE) ## hapster
# cnv.fn <- list.files("/mctp/users/noshadh/work/HLACNV/data/cohort_177/cohort_177_cnvex", full.names = TRUE) ## cnvex
#
# mctp_info_normal.all[, calls := paste0(tstrsplit(fasta, "refs/")[[1]], "calls")]
# mctp_info_normal.all = readRDS("/mctp/users/noshadh/data/polytec/mctp-covs-normal-updated-with-cnv.rds")
# mctp_info_normal = copy(mctp_info_normal.all)
# # mctp_info_normal = mctp_info_normal[capture == "Onco1500-v4"]
# mctp_info_normal = mctp_info_normal[capture == "Agilent-v4"]
# agi.digest = mctp_info_normal[!is.na(mctp_info_normal$digest)]
# #
# cnv.fn = mctp_info_normal$cnv
# cnv.fn <- sort(cnv.fn)
# hapster.fn = mctp_info_normal$calls
# hapster.fn <- unlist(strsplit(hapster.fn, "/calls"))
# hapster.fn <- sort(hapster.fn)
# hapster.fn <- hapster.fn[(which(table(hapster.fn) <2))]
#
# cnv.fn     <- cnv.fn[1:200]
# hapster.fn <- hapster.fn[1:200]
# HLACNV.opts <- MHCnvex::.HLACNV.getOpts()
# HLACNV.opts$haplotyper.mode = "polytect"
# samp = 10
# fasta = mctp_info_normal$fasta[samp]
# calls = mctp_info_normal$calls[samp]
# bam.n = mctp_info_normal$bam.n[samp]
# cnv <- mctp_info_normal$cnv[samp]
# pool <- HLACNV.makePool(hapster.fn, cnv.fn, 20, HLACNV.opts)
#

# i=20
# mctp_info_normal.all = readRDS("/mctp/users/noshadh/data/polytec/mctp-covs-normal-updated-with-cnv.rds")
# mctp_info_normal.all[, calls := paste0(tstrsplit(fasta, "refs/")[[1]], "calls")]
# mctp_info_normal = copy(mctp_info_normal.all)
# # mctp_info_normal = mctp_info_normal[capture == "Onco1500-v4"]
# mctp_info_normal = mctp_info_normal[capture == "Agilent-v4"]
# agi.digest = mctp_info_normal[!is.na(mctp_info_normal$digest)]
#
#
#
# i = which(agi.digest$patient == "TP_2114")
# digest = readRDS(agi.digest$digest[i])
# cnv = readRDS(agi.digest$cnv[i])
# sample_id = agi.digest$patient[i]
# n.bam = agi.digest$bam.n[i]
# t.bam = agi.digest$bam.t[i]
# fasta = agi.digest$fasta[i]
# calls = agi.digest$calls[i]
# purity = digest$model.purity
# ploidy = digest$model.ploidy
# digest = agi.digest$digest[i]
# cnv = agi.digest$cnv[i]

# sample_id = "MO_1042"
# n.bam = "/mctp/users/noshadh/data/polytec/MCTP_vcfs/MO_1042/alignments/SI_5028-C0LARACXX-7-All-merged.bam"
# t.bam = "/mctp/users/noshadh/data/polytec/MCTP_vcfs/MO_1042/alignments/SI_5029-C0LARACXX-7_all_germline_imputed.bam"
# fasta = "/mctp/users/noshadh/data/polytec/MCTP_vcfs/MO_1042/refs/MO_1042_all_germline_imputed.fa"
# calls = "/mctp/users/noshadh/data/polytec/MCTP_vcfs/MO_1042/calls"
# cnv   = "/mctp/share/mioncoseq2/cnv/MO_1042-SI_5029.541-SI_5028.531/MO_1042-SI_5029.541-SI_5028.531.rds"
#


# test for current state of the code:
# i=20
# mctp_info_normal.all = readRDS("/mctp/users/noshadh/data/polytec/mctp-covs-normal-updated-with-cnv.rds")
# mctp_info_normal.all[, calls := paste0(tstrsplit(fasta, "refs/")[[1]], "calls")]
# mctp_info_normal = copy(mctp_info_normal.all)
# # mctp_info_normal = mctp_info_normal[capture == "Onco1500-v4"]
# mctp_info_normal = mctp_info_normal[capture == "Agilent-v4"]
# agi.digest = mctp_info_normal[!is.na(mctp_info_normal$digest)]
#
#
#
# i = which(agi.digest$patient == "TP_2114")
# digest = readRDS(agi.digest$digest[i])
# cnv = readRDS(agi.digest$cnv[i])
# sample_id = agi.digest$patient[i]
# n.bam = agi.digest$bam.n[i]
# t.bam = agi.digest$bam.t[i]
# fasta = agi.digest$fasta[i]
# calls = agi.digest$calls[i]
# purity = digest$model.purity
# ploidy = digest$model.ploidy
# digest = agi.digest$digest[i]
# cnv = agi.digest$cnv[i]



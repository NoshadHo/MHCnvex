.HLACNV.addPeak <- function(exons, HLACNV.opts, peak.alignid=NULL) {
    ## if peak.alignid availabe
    if (!is.null(peak.alignid)) {
      col.order        <- colnames(exons)
      exons            <- merge.data.table(exons, peak.alignid, by = c("align.id", "hla.type"), all.x = TRUE, sort = FALSE)
      exons$peak.PASS  <- ifelse(is.na(exons$peak.id), FALSE, TRUE)
      exons$cov.PASS   <- ifelse(is.na(exons$peak.id), FALSE, TRUE) ## if selected from pool, it's coverage shuold be sufficient
      ## TODO: make sure it makes sense to change the cov.pass
      return(exons)
    }
    
    if(any(!is.na(exons$BAF.n))) { ## if normal sample available
      cov.tbl <- exons[, .(cov=hap1.n.cov+hap2.n.cov, start=start, hla.type=hla.type, cov.PASS=cov.PASS)]
    } else {
      ## assumption: each gene cnv affect the whole gene similarly
      cov.tbl <- exons[, .(cov=hap1.t.cov+hap2.t.cov, start=start, hla.type=hla.type, cov.PASS=cov.PASS)]
    }
    
    ## finding peaks
    peaks <- .HLACNV.detectPeaks(exons, cov.tbl, HLACNV.opts, w=1)
    
    ## filtering
    peaks <- .HLACNV.filterPeaks(peaks, cov.tbl, HLACNV.opts)
    
    exons$peak.PASS  <- peaks$peak.adj.PASS
    exons$peak.id  <- peaks$peak.id
    
    return(exons)
}


.HLACNV.detectPeaks <- function(exons, cov.tbl, HLACNV.opts, w=1) { # source: http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
    cov.tbl.split <- split(cov.tbl, cov.tbl$hla.type)
    peaks  <- lapply(cov.tbl.split, function(gene.covs) {

        peaks      <- .HLACNV.detectPeaksInner(x=gene.covs$start, y=gene.covs$cov, w=1, HLACNV.opts)
        peaks.tbl  <- data.table(gene = gene.covs$hla.type[1], peak.loc = peaks$x, peak.index = peaks$i,
                                 peaks.y.smooth = peaks$y.hat[peaks$i],
                                 peak.cov = gene.covs[peaks$i]$cov, peak.cov.PASS = gene.covs[peaks$i]$cov.PASS,
                                 id=paste(gene.covs$hla.type[1],peaks$x, sep="-"))
        return(peaks.tbl)
    })
    peaks <- rbindlist(peaks)

    return(peaks)
}

.HLACNV.detectPeaksInner <- function(x, y, w=1, HLACNV.opts) { # perform smoothing for each sample as well
    ## smooth y
    y.smooth <- loess(y ~ x, span = HLACNV.opts$filter.perBase.loess.smooth)$fitted

    ## detecting
    n <- length(y.smooth)
    y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
    i.max <- which(delta <= 0) + w
    list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

.HLACNV.basepair2Peak <- function(exons, HLACNV.opts, align.basepair=NULL) {
    ## find highest coverage peak of each gene and handle regions with no peak
    missed.peaks <- exons[, .(
        start=min(start), end = max(end),
        hap1.n.cov.raw = mean(hap1.n.cov.raw), hap2.n.cov.raw = mean(hap2.n.cov.raw),
        hap1.t.cov.raw = mean(hap1.t.cov.raw), hap2.t.cov.raw = mean(hap2.t.cov.raw),
        hap1.n.cov = mean(hap1.n.cov), hap2.n.cov = mean(hap2.n.cov),
        hap1.t.cov = mean(hap1.t.cov), hap2.t.cov = mean(hap2.t.cov),
        hetero     = hetero[1]

    ), by = .(hla.type, hap1,hap2, peak.id)]
    missed.peaks <- missed.peaks[, .SD[hap1.n.cov+hap2.n.cov == max(hap1.n.cov+hap2.n.cov)], by = hla.type]
    missed.peaks <- missed.peaks[is.na(peak.id)]
    missed.peaks[is.na(peak.id), peak.id := sprintf("%s-7410CCC", hla.type)]

    ## aggregating each peak as one tile
    peakTile <- exons[peak.PASS == TRUE & cov.PASS == TRUE, .(
        start=min(start), end = max(end),
        hap1.n.cov.raw = mean(hap1.n.cov.raw), hap2.n.cov.raw = mean(hap2.n.cov.raw),
        hap1.t.cov.raw = mean(hap1.t.cov.raw), hap2.t.cov.raw = mean(hap2.t.cov.raw),
        hap1.n.cov = mean(hap1.n.cov), hap2.n.cov = mean(hap2.n.cov),
        hap1.t.cov = mean(hap1.t.cov), hap2.t.cov = mean(hap2.t.cov),
        hetero     = hetero[1]

    ), by = .(peak.id,hap1,hap2,hla.type)]

    ## if any gene is missing a peak, we shuold add one
    peakTile <- rbind(peakTile, missed.peaks)

    ## calc new lbaf for peaktiles (and dismissing basewise baf)
    peakTile[, BAF.n := hap1.n.cov/(hap1.n.cov+hap2.n.cov+0.001)]
    peakTile[, BAF.t := hap1.t.cov/(hap1.t.cov+hap2.t.cov+0.001)]


    ## filters
    if(any(!is.na(peakTile$hap1.n.cov) | !is.na(peakTile$hap2.n.cov))) { ## if normal sample available
        peakTile[, baf.PASS  := abs(BAF.n - 0.5) < HLACNV.opts$filter.perBase.max.n.dev]
    } else {
        peakTile[, baf.PASS  := TRUE]
    }

    ## used clean basepair, so filters are passed
    peakTile[, cov.PASS  := TRUE]
    peakTile[, peak.PASS := TRUE]

    return(peakTile)
}


.HLACNV.addPeak2BasePair <- function(basepair, HLACNV.opts, exons=NULL) {
  ## if no exons, or exon is in peakTile mode (peakBase is required)
  if (is.null(exons) | (all(exons$peak.PASS == TRUE) & all(table(exons$peak.id) == 1))) {
    exons <- .HLACNV.makeExonsPerBase(basepair, ref, "peakBase", HLACNV.opts)
  }
  ## given the logic of peak finding, calculating exon is essential (needs total coverage at each point)
  
  ## add peak
  basepair.split <- split(basepair, basepair$seqnames)
  subbps <- lapply(basepair.split, function(subbp) {
    subexons <- exons[hla.type == subbp$hla.type[1]]
    hap <- ifelse(subbp$seqnames[1] %in% subexons$hap1[1], "hap1", "hap2")
    subbp    <- merge.data.table(subbp, subexons[, c(..hap, "align.id", "peak.PASS", "peak.id", "hetero")], 
                                 by.x = c("seqnames", "align.id"), by.y = c(hap, "align.id"), all.x = TRUE)
    subbp[is.na(peak.PASS), peak.PASS := FALSE]
  })
  basepair.new <- rbindlist(subbps)
  
  ## estimate the homozygous alleles
  basepair.new.split <- split(basepair.new, basepair$hla.type)
  subbps <- lapply(basepair.new.split, function(subbp) {
    if(all(subbp$hetero)) {
      subbp[, allele.type := "heterozygous"]
    } else {
      alleles     <- subbp[peak.PASS == TRUE, .(cov = mean(n.cov)), by = seqnames]
      maj.allele  <- alleles[which.max(alleles$cov)]$seqnames
      min.allele  <- alleles[which.max(alleles$cov)]$seqnames
      
      subbp[seqnames == maj.allele, allele.type := "homozygous-maj"]
      subbp[seqnames == min.allele, allele.type := "homozygous-min"]
    }
    
    return(subbp)
  })
  
  basepair.final <- rbindlist(subbps)
  
  return(basepair.final)
}
# 
# .HLACNV.exon2peakCoverage <- function(peakTile, ref) {
#   ## allele1 imp coordinate
#   allele1.rawCoords <- merge.data.table(peakTile[, c("peak.id", "hap1", "hap2", "hla.type", "start", "end")], ref$mapping.tbl[, c("seqnames", "imp.id", "align.id")],
#                                         by.x = c("hap1", "start"), by.y = c("seqnames", "align.id"))
#   setnames(allele1.rawCoords, "imp.id", "start.imp.id")
#   allele1.rawCoords <- merge.data.table(allele1.rawCoords, ref$mapping.tbl[, c("seqnames", "imp.id", "align.id")],
#                                         by.x = c("hap1", "end"), by.y = c("seqnames", "align.id"))
#   setnames(allele1.rawCoords, "imp.id", "end.imp.id")
#   
#   ## allele2 imp coordinate
#   allele2.rawCoords <- merge.data.table(peakTile[, c("peak.id", "hap1", "hap2", "hla.type", "start", "end")], ref$mapping.tbl[, c("seqnames", "imp.id", "align.id")],
#                                         by.x = c("hap2", "start"), by.y = c("seqnames", "align.id"))
#   setnames(allele2.rawCoords, "imp.id", "start.imp.id")
#   allele2.rawCoords <- merge.data.table(allele2.rawCoords, ref$mapping.tbl[, c("seqnames", "imp.id", "align.id")],
#                                         by.x = c("hap2", "end"), by.y = c("seqnames", "align.id"))
#   setnames(allele2.rawCoords, "imp.id", "end.imp.id")
#   
#   ## calculate the coverage
#   rawcoords <- data.table(seqnames = c(allele1.rawCoords$hap1, allele2.rawCoords$hap2),
#                           start    = c(allele1.rawCoords$start.imp.id, allele2.rawCoords$start.imp.id),
#                           end      = c(allele1.rawCoords$end.imp.id, allele2.rawCoords$end.imp.id),
#                           peak.id  = c(allele1.rawCoords$peak.id, allele2.rawCoords$peak.id))
#   allele1.exon <- MHCnvex:::.HLACNV.importBam(ref$t.bam, ref$n.bam, rawcoords, ref$fasta, HLACNV.opts)
#   allele1.exon <- MHCnvex:::.HLACNV.normCoverageBasepair()
#   allele1.exon.split  <- split(allele1.exon, allele1.exon$seqnames)
#   
#   allele1.exon.split <- 
#   
#   
# }
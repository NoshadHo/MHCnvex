.HLACNV.normGCLr <- function(cnv, hla.exons, HLACNV.opts, opts) { ##  remove opts from the options later on and add required options to hlacnv.opts
    exons.obj <- copy(hla.exons)
    cnv.hla   <- copy(cnv)
    cnv.tile  <- as.data.table(cnv.hla$tile)
    hla.tile  <- as.data.table(exons.obj)
    hla.tile[, ':='(n.cov = hap1.n.cov+hap2.n.cov, t.cov=hap1.t.cov+hap2.t.cov)]
    
    tile.dt   <- GRanges(seqnames = c(cnv.tile$seqnames, hla.tile$hla.type, hla.tile$hap1, hla.tile$hap2), 
                         ranges = IRanges(c(cnv.tile$start, rep(1, nrow(hla.tile)), rep(1, nrow(hla.tile)), rep(1, nrow(hla.tile))), 
                                          c(cnv.tile$end, rep(1, nrow(hla.tile)), rep(1, nrow(hla.tile)), rep(1, nrow(hla.tile)))),
                         n.cov = c(cnv.tile$n.cov, hla.tile$n.cov, hla.tile$hap1.n.cov, hla.tile$hap2.n.cov), 
                         t.cov = c(cnv.tile$t.cov, hla.tile$t.cov, hla.tile$hap1.t.cov, hla.tile$hap2.t.cov),
                         gc = c(cnv.tile$gc, hla.tile$gc, hla.tile$gc, hla.tile$gc), 
                         hq = c(cnv.tile$hq, rep(TRUE, nrow(hla.tile)), rep(TRUE, nrow(hla.tile)), rep(TRUE, nrow(hla.tile))), 
                         target = c(cnv.tile$target, rep(TRUE, nrow(hla.tile)), rep(TRUE, nrow(hla.tile)), rep(TRUE, nrow(hla.tile))),
                         type = c(rep("cnv", nrow(cnv.tile)), rep("hla", nrow(exons.obj)), rep("hap1", nrow(exons.obj)), rep("hap2", nrow(exons.obj))))
    cnv.hla$tile <- tile.dt
    
    cnv.hla <- suppressWarnings(cnvex::addLogRatio(cnv.hla, NULL, opts))
    
    exons.obj$lr.dip <- cnv.hla$tile$tn.lr[cnv.hla$tile$type == "hla"]
    exons.obj$lr.hap1 <- cnv.hla$tile$tn.lr[cnv.hla$tile$type == "hap1"]
    exons.obj$lr.hap2 <- cnv.hla$tile$tn.lr[cnv.hla$tile$type == "hap2"]
    cnv.hla$tile     <- cnv.hla$tile[cnv.hla$tile$type == "cnv"]
    
    return(list(exons.obj, cnv.hla))
}
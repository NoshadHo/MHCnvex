.HLACNV.makeMismatchPerBase <- function(basepair, ref, HLACNV.opts, peak.alignid = NULL) {   ## TODO: the peaktile logic is problematic, fix it
  exons <- .HLACNV.makeExonsPerBase(basepair, ref, "peakBase", HLACNV.opts, peak.alignid)
  peak.loci <- exons[!is.na(peak.id)]
  peak.loci[, mismatch := seq1 != seq2]
  mismatch <- peak.loci[mismatch == TRUE]
  mismatch$group <- paste0(mismatch$hla.type, "-HLA")
  return(mismatch)
  
}

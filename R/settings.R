.HLACNV.getConsts <- function() {
    ## to generate this HLA positions, if any of the HLAs had overlap with the other one,
    # I manually changes the second gene start point to be start(next) = end(previous)+2

    HLA.CORDS <- system.file("extdata", "GRCh38/GRCh38.HLA.pos4.rds", package = "MHCnvex")
    return(readRDS(HLA.CORDS))
}

.HLACNV.getExternalRefs <- function(ref.name) {
    # ref.name: exon.tiles, align.seqs.polytect, align.seqs.hapster
    if (ref.name == "polytect.exons") {
        exon.tiles <- system.file("extdata", "GRCh38/polytect_alt_genes.bed", package = "MHCnvex")
        external.ref <- data.table::fread(exon.tiles, col.names = c("seqnames", "start", "end", "idx"), select = c(1,2,3,7))
    } else if(ref.name == "hapster.exons") {
        exon.tiles <- system.file("extdata", "GRCh38/hapster_alt_genes.bed", package = "MHCnvex")
        external.ref <- data.table::fread(exon.tiles, col.names = c("seqnames", "start", "end", "idx"), select = c(1,2,3,7))
    } else if(ref.name == "align.seqs.polytect") {
        align.seqs.polytect <- system.file("extdata", "haplotyper/all_seqs_polytect.txt", package = "MHCnvex")     # aligned seqs
        external.ref <- data.table::fread(align.seqs.polytect,sep = "\t",col.names = c("names","seqs"),header = FALSE)
    } else if(ref.name == "align.seqs.hapster") {
        align.seqs.hapster <- system.file("extdata", "haplotyper/all_seqs_hapster.txt", package = "MHCnvex")     # aligned seqs
        external.ref <- data.table::fread(align.seqs.hapster,sep = "\t",col.names = c("names","seqs"),header = FALSE)
    }

    return(external.ref)
}

.HLACNV.getOpts <- function() {
    HLACNV.opts <- system.file("extdata", "settings/MHCopts.R", package = "MHCnvex")
    ENV <- new.env(parent = .BaseNamespaceEnv)
    source(HLACNV.opts, local=ENV)
    HLACNV.opts <- ENV$HLACNV.opts
    return(HLACNV.opts)
}

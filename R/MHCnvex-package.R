#' MHCnvex
#'
#' @name MHCnvex
#' @docType package
#'
#' @section main.R:
#' @importFrom data.table fwrite
#' @importFrom ggplot2 ggsave ggtitle
#' @importFrom GenomicRanges width
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @section blend.R:
#' @importFrom GenomicRanges findOverlaps mcols makeGRangesFromDataFrame
#' @importFrom data.table copy data.table fcase
#'
#' @section exon.R:
#' @importFrom data.table copy data.table merge.data.table rbindlist tstrsplit
#' @importFrom Biostrings DNAString letterFrequency
#'
#' @section filter.R:
#' @importFrom stats sd
#'
#' @section gc.R:
#' @importFrom IRanges IRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table copy as.data.table
#'
#' @section mismatch.R:
#' @importFrom data.table copy merge.data.table
#'
#' @section models.R:
#' @importFrom GenomicRanges seqnames end start sort findOverlaps makeGRangesFromDataFrame mcols<-
#' @importFrom data.table as.data.table foverlaps merge.data.table rbindlist setkey rleid
#'
#' @section output.R:
#' @importFrom stats median
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom data.table as.data.table data.table merge.data.table rbindlist
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @section peak.R:
#' @importFrom data.table data.table foverlaps rbindlist setkey
#' @importFrom stats loess
#' @importFrom zoo rollapply zoo
#'
#' @section plots.R:
#' @importFrom GenomicRanges end start
#' @importFrom data.table as.data.table data.table tstrsplit
#' @importFrom ggplot2 aes element_blank geom_hline geom_label geom_point ggplot scale_color_manual scale_y_continuous theme theme_void scale_size_manual geom_segment
#' @importFrom ggpubr theme_pubr
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork inset_element
#' @importFrom stringr str_pad
#' @importFrom matrixStats colMedians
#' @importFrom ComplexHeatmap Heatmap decorate_heatmap_body HeatmapAnnotation draw
#' @importFrom grid grid.lines gpar
#'
#' @section ref.R:
#' @importFrom GenomicRanges seqnames width seqnames<-
#' @importFrom GenomeInfoDb genome genome<- seqinfo
#' @importFrom DelayedArray rowRanges
#' @importFrom VariantAnnotation fixed geno readVcf ScanVcfParam
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table data.table as.data.table fread merge.data.table rbindlist foverlaps
#'
#' @section tile.R:
#' @importFrom BiocGenerics get
#' @importFrom Biostrings letterFrequency getSeq
#' @importFrom rtracklayer export
#' @importFrom data.table fread rbindlist
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges
#' @importFrom GenomeInfoDb sortSeqlevels
#'
#' @section tilebaf.R:
#' @importFrom stats na.omit
#'
#'
NULL

# NCmisc::list.functions.in.file("./R/tile.R")

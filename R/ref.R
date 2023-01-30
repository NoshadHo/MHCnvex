.HLACNV.readRef <- function(fasta, calls ,HLACNV.opts) {
    # alignrefs.fn: folder containing all the refs
    # sequence
    ref.stringset <- readDNAStringSet(fasta)
    ref.imputed    <- data.table(name = names(ref.stringset), seq = as.character(ref.stringset), width = width(ref.stringset),
                                 hla = do.call(rbind,strsplit(names(ref.stringset), "\\*"))[,1])
    ref.align     <- .HLACNV.readAlignSeqs(fasta, HLACNV.opts)
    ref.raw       <- .HLACNV.makeSeqRaw(ref.align)
    # seqInfo
    ref.info <- seqinfo(ref.stringset)
    genome(ref.info) = "hg38"

    ## alignment
    mapping.tbl <- makeMappingTbl(ref.imputed, ref.raw, ref.align, calls, HLACNV.opts)
    mapping.tbl <- .HLACNV.mapExon2Ref(mapping.tbl)
    ## exon regions
    ## dev state, needs to be teste
    # tile.raw <- .HLACNV.makeTileRaw(ref.imputed, HLACNV.refs)
    # ref.tbl  <- .HLACNV.tileRaw2Imp(tile.raw, mapping.tbl)
    # tile.imp <- .HLACNV.annotatetile(ref.tbl$imp.exons, ref.stringset, HLACNV.opts)


    ref <- list(seq=ref.stringset, seqi=ref.info, seq.imputed=ref.imputed, seq.align=ref.align, seq.raw = ref.raw, mapping.tbl = mapping.tbl,
                fasta=fasta)


    return(ref)
}


.HLACNV.readAlignSeqs <- function(fasta, HLACNV.opts) {
    if (HLACNV.opts$haplotyper.mode == "polytect") {
        seqs <- .HLACNV.getExternalRefs("align.seqs.polytect")
    } else if(HLACNV.opts$haplotyper.mode == "hapster") {
        seqs <- .HLACNV.getExternalRefs("align.seqs.hapster")
    }

    # seqs           <- fread(align.seqs.fn,sep = "\t",col.names = c("names","seqs"),header = FALSE)
    ref.stringset  <- readDNAStringSet(fasta)
    seqs.fasta     <- seqs[names %in% names(ref.stringset)]
    seqs.fasta$hla <- do.call(rbind,strsplit(names(ref.stringset), "\\*"))[,1]
    return(seqs.fasta)
}


.HLACNV.makeSeqRaw <- function(ref.align) {
    seq.raw <- lapply(1:nrow(ref.align), function(i) {
        row <- ref.align[i,]
        seq <- strsplit(row$seqs, "")[[1]]
        seq <- seq[seq %in% c("A","T","C","G")]
        row$seq <- paste0(seq, collapse = "")
        return(row)
    })
    seq.raw <- rbindlist(seq.raw)
    return(seq.raw)
}

makeMappingTbl <- function(ref.imputed, ref.raw, ref.align, calls, HLACNV.opts) {

    mapping.Raw2Align <- .HLACNV.mapRaw2Align(ref.raw, ref.align)
    mapping.Imp2Raw   <- .HLACNV.mapRaw2Imputed(ref.imputed, calls, HLACNV.opts)
    mapping.tbl <- merge.data.table(mapping.Imp2Raw, mapping.Raw2Align[,c("seq","raw.id","align.id", "seqnames")], by=c("seqnames","raw.id"),all = FALSE)
    if(sum(mapping.tbl$seq.x != mapping.tbl$seq.y, na.rm=T) == 0) {
        print("Mapping filter: PASS!")
    }
    mapping.tbl <- keepAlignedMaps(mapping.tbl)
    colnames(mapping.tbl)[3] = "seq.raw"
    return(mapping.tbl[,c("seqnames", "imp.id", "raw.id", "align.id","seq.raw", "hla.type")])

}

keepAlignedMaps <- function(mapping.tbl) {
    ## any base not aligned, goes straight to ... :)
    mapping.tbl.split <- split(mapping.tbl, mapping.tbl$hla.type)
    mapping.tbl.inner <- lapply(mapping.tbl.split, function(subtable) {
        subtable    <- split(subtable, subtable$seqnames)
        common.pos  <- intersect(subtable[[1]]$align.id, subtable[[2]]$align.id)
        subtable.in <- rbind(subtable[[1]][align.id %in% common.pos], subtable[[2]][align.id %in% common.pos])
    })
    return(rbindlist(mapping.tbl.inner))
}

#
.HLACNV.mapRaw2Align <- function(seq.raw, seq.align) { ##############

    mapping <- lapply(1:nrow(seq.raw), function(i) {
        # print(seq.raw[i,]$names)
        seq.raw.sel         <- data.table(seq = (strsplit(seq.raw[i,]$seq, ""))[[1]], raw.id = 1:nchar(seq.raw[i,]$seq))
        seq.align.sel       <- data.table(seq = (strsplit(seq.align[i,]$seq, ""))[[1]], align.id = 1:nchar(seq.align[i,]$seq))

        seq.align.sel        <- seq.align.sel[seq %in% c("T","C","G","A")]
        seq.align.sel$raw.id <- 1:nrow(seq.align.sel)

        seq.align.sel$seqnames <- seq.raw[i,]$name
        seq.align.sel$hla.type <- seq.raw[i,]$hla

        return(seq.align.sel)
    })
    mapping <- rbindlist(mapping)
    return(mapping) }


.HLACNV.mapRaw2Imputed <- function(ref.imputed, calls, HLACNV.opts) {
    if (HLACNV.opts$haplotyper.mode == "polytect") {
        ref.alts <- .HLACNV.imputeGermlineAltsPolytect(ref.imputed, calls, HLACNV.opts)
    } else if(HLACNV.opts$haplotyper.mode == "hapster") {
        ref.alts <- .HLACNV.imputeGermlineAltsHapster(ref.imputed, calls, HLACNV.opts)
    }


    mapping.table <- lapply(1:nrow(ref.imputed), function(i) {
        row          <- ref.imputed[i,]
        alts         <- as.data.table(ref.alts[as.character(seqnames(ref.alts)) %in% row$name])
        seq          <- data.table(seq = strsplit(row$seq,"")[[1]], imp.id = NA_real_)
        seq$imp.id   <- 1:nrow(seq)
        if (nrow(alts) == 0) { # don't fix anything if no alterations
            seq[, raw.id := 1:.N]
            seq$seqnames <- row$name
            seq$hla.type <- row$hla
            return(seq)
        }
        for(j in 1:nrow(alts)) { ## replacing alt with ref
            var <- alts[j]
            alt <- var$ALT
            ref <- var$REF
            # str_sub(seq, var$start, var$start+alt.len-1) <- var$REF

            ##Looks fancy, but just ID part is just to make sure the order of ids would be correct
            ref.dt <- data.table(seq = strsplit(ref,"")[[1]], imp.id = (seq[(var$start-1)]$imp.id+0.01*(1:length(strsplit(ref,"")[[1]]))))

            seq <- rbind(seq[1:(var$start-1)], ref.dt, seq[(var$start+nchar(alt)):(.N),])
        }
        seq[, raw.id := 1:.N]
        seq$seqnames <- row$name
        seq$hla.type  <- row$hla
        return(seq)
    })
    mapping.table <- rbindlist(mapping.table)
    return(mapping.table)
}


## old version
.HLACNV.imputeGermlineAltsPolytect <- function(ref.imputed, calls, HLACNV.opts) {
    param <-ScanVcfParam(
        info=c("AF"),
        geno=c("AD")
    )
    ref.alts <- lapply(unique(ref.imputed$hla), function(hla) {

        # row       <- ref.imputed[i,]
        call      <- list.files(path = sprintf("%s/",calls), full.names = TRUE, pattern = sprintf("*_%s_germline_filtered.vcf", hla))
        vcf       <- suppressWarnings(readVcf(call,param = param))
        ## if more than 1 alt, polytect won't impute, so doesn't matter
        vcf       <- vcf[which(lengths(fixed(vcf)$ALT) == 1)]
        var.g     <- rowRanges(vcf)
        seqlevels(var.g) <- unique(as.character(seqnames(var.g)))
        
        var.g$REF <- as.character(var.g$REF)
        # var.g$ALT <- as.character(unlist(var.g$ALT))
        var.g$ALT <- sapply(var.g$ALT, function(alt) {
            alt <- alt[1]
            return(as.character(unlist(alt)))

        })
        AD        <- do.call(rbind, geno(vcf)$AD[,1])[,1:2]
        if (is.null(AD)) {
            var.g$AD <- matrix(numeric(),ncol=2)
        } else {
            var.g$AD <- matrix(AD, ncol = 2,)
        }
        # var.g$AD  <- t(sapply(geno(vcf)$AD, "[", 1:2))
        var.g$ADF <- var.g$AD[,2]/(var.g$AD[,1] + var.g$AD[,2]) ## alt read support
        var.g$ADF <- ifelse(is.na(var.g$ADF), 0, var.g$ADF)
        print(paste0(hla, ": ",length(var.g[var.g$ADF > HLACNV.opts$ref.gemline.impute.thr])))
        return(var.g)
    })
    ref.alts <- suppressWarnings(do.call(c, ref.alts))
    ref.alts <- ref.alts[ref.alts$ADF > HLACNV.opts$ref.gemline.impute.thr]
    return(ref.alts)
}

## newer version, but needs to be modified to match old version
.HLACNV.imputeGermlineAltsHapster <- function(ref.imputed, calls, HLACNV.opts) {
    param <-ScanVcfParam(
        info=c("AF"),
        geno=c("AD")
    )

    call      <- list.files(path = sprintf("%s/",calls), full.names = TRUE, pattern = "*_germline_filtered.vcf$")
    vcf       <- suppressWarnings(readVcf(call,param = param))
    ## if more than 1 alt, Hapster won't impute, so doesn't matter
    vcf       <- vcf[which(lengths(fixed(vcf)$ALT) == 1)]
    var.g     <- rowRanges(vcf)
    var.g$REF <- as.character(var.g$REF)
    # var.g$ALT <- as.character(unlist(var.g$ALT))
    var.g$ALT <- sapply(var.g$ALT, function(alt) {
        alt <- alt[1]
        return(as.character(unlist(alt)))

    })
    AD        <- do.call(rbind, geno(vcf)$AD[,1])[,1:2]
    if (is.null(AD)) {
        var.g$AD <- matrix(numeric(),ncol=2)
    } else {
        var.g$AD <- matrix(AD, ncol = 2,)
    }
    # var.g$AD  <- t(sapply(geno(vcf)$AD, "[", 1:2))
    var.g$ADF <- var.g$AD[,2]/(var.g$AD[,1] + var.g$AD[,2]) ## alt read support
    var.g$ADF <- ifelse(is.na(var.g$ADF), 0, var.g$ADF)
    print(paste0("Gemline filter: ",length(var.g[var.g$ADF > HLACNV.opts$ref.gemline.impute.thr])))

    var.g <- var.g[var.g$ADF > HLACNV.opts$ref.gemline.impute.thr]
    return(var.g)
}

.HLACNV.mapExon2Ref <- function(mapping.tbl) {
    if (HLACNV.opts$haplotyper.mode == "polytect") {
        exons <- .HLACNV.getExternalRefs("polytect.exons")
    } else if(HLACNV.opts$haplotyper.mode == "hapster") {
        exons <- .HLACNV.getExternalRefs("hapster.exons")
    }

    exons.sample <- exons[seqnames %in% unique(mapping.tbl$seqnames)]
    setkey(exons.sample, seqnames, start, end)
    hits <- foverlaps(mapping.tbl[, .(seqnames = seqnames ,start = raw.id, end = raw.id)],
                      exons.sample, by.x = c("seqnames", "start", "end"), which = TRUE)
    map.id <- cbind(mapping.tbl[hits$xid],
                    exons.sample[hits$yid, .(is.exon = ifelse(is.na(hits$yid), FALSE, TRUE), exon.start = start, exon.end = end)])

    map.id[is.exon == TRUE, ref.exon.id := (paste0(hla.type,"-",exon.start)), by = .(seqnames, exon.start, exon.end)]

    return(map.id)
}

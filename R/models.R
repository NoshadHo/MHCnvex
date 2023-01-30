.HLACNV.segmentation <- function(mhla, opts) {
    # opts <- cnvex::getOpts(system.file(sprintf("extdata/settings/%s.R", "exome-pair"), package="cnvex"))
    tile <- mhla$tile[seqnames(mhla$tile) == "chr6"]
    # tile <- mhla$tile
    tile$arm <- factor(tile$arm, levels = unique(tile$arm))
    seg  <- cnvex:::jointSegment(tile, opts)
    opts$prune.nvar <- FALSE
    seg  <- cnvex:::pruneSegments(seg, tile, mhla$stats, opts)
    return(seg)
}

.HLACNV.modelHLA <- function(exons, mismatch, tile, cnv, HLACNV.opts, opts) {
    mcnv.dummy  <- cnvex:::modelCnv("tumor", cnv, NULL, NULL, opts) ## use stats and sex information of mcnv.dummy
    ####### VAR
    # var.raw   <- .HLACNV.blend.cnvVar.hla(mismatch, cnv, HLACNV.opts)
    chr6.vars   <- as.data.table(.HLACNV.blend.cnvVar.hla(mismatch, cnv, HLACNV.opts))
    # chr6.coords <- c(start = min(var.raw[seqnames(var.raw) == "chr6"]$idx),
                     # end   = max(var.raw[seqnames(var.raw) == "chr6"]$idx))
    # chr6.vars  <- as.data.table(var.raw[var.raw$idx > chr6.coords['start'] & var.raw$idx < chr6.coords['end']])
    chr6.vars$seqnames <- as.character(chr6.vars$seqnames)
    ## mark hla genes
    chr6.vars[group != "GRCh38"]$group <- sprintf("%s-HLA",as.character((chr6.vars[chr6.vars$group != "GRCh38"]$seqnames)))
    chr6.vars$AF <- chr6.vars$t.AF
    chr6.vars$DP <- chr6.vars$t.DP

    ###### TILE
    ## seperate chr6
    # tile.raw    <- .HLACNV.blend.cnvTile.hla(exons, cnv, HLACNV.opts)
    # tile.raw$arm<- as.character(tile.raw$arm)
    chr6.tiles     <- as.data.table(.HLACNV.blend.cnvTile.hla(exons, cnv, HLACNV.opts))
    chr6.tiles$arm <- as.character(chr6.tiles$arm)
    # chr6.coords <- c(start = min(tile.raw[seqnames(tile.raw) == "chr6"]$idx),
                     # end   = max(tile.raw[seqnames(tile.raw) == "chr6"]$idx))
    # chr6.tiles  <- as.data.table(tile.raw[tile.raw$idx > chr6.coords['start'] & tile.raw$idx < chr6.coords['end']])
    chr6.tiles$seqnames <- as.character(chr6.tiles$seqnames)
    ## mark hla genes
    chr6.tiles[group != "GRCh38"]$group <- sprintf("%s-HLA",as.character((chr6.tiles[chr6.tiles$group != "GRCh38"]$seqnames)))
    chr6.tiles <- .HLACNV.add.baf(chr6.tiles[target == TRUE & hq == TRUE],
                                  chr6.vars[!(mask.loose | mask.strict)==TRUE & PASS==TRUE], opts)



    ## fix seqnames and positions
    syncdata <- .HLACNV.syncTileVar(chr6.tiles, chr6.vars, HLACNV.opts)

    ## for wgs
    if(HLACNV.opts$panel.mode == "wgs") {
        syncdata <- .HLACNV.makeWgsCompatible(syncdata, HLACNV.opts)
    }

    mcnv.dummy$tile     <- makeGRangesFromDataFrame(syncdata[['tile']], keep.extra.columns = TRUE)
    mcnv.dummy$var      <- makeGRangesFromDataFrame(syncdata[['var']], keep.extra.columns = TRUE)

    ## adding hap names
    mcnv.dummy$haps <- exons[,.(hap1 = hap1[1], hap2=hap2[1], group = paste0(hla.type[1], "-HLA"), hetero = hetero[1]),by = hla.type]


    ## filtering
    mcnv.dummy$tile <- mcnv.dummy$tile[mcnv.dummy$tile$target == TRUE & mcnv.dummy$tile$hq == TRUE]
    mcnv.dummy$var  <- mcnv.dummy$var[!(mcnv.dummy$var$mask.loose | mcnv.dummy$var$mask.strict) & mcnv.dummy$var$PASS]

    mcnv.dummy$tile <- GenomicRanges::sort(mcnv.dummy$tile)
    mcnv.dummy$var  <- GenomicRanges::sort(mcnv.dummy$var)

    ## sort arms
    # mcnv.dummy$tile$arm <- factor(mcnv.dummy$tile$arm, levels = levels(cnv$tile$arm), ordered = TRUE)

    return(mcnv.dummy)

}

## ~3Mbp distance between HLA_A and HLA_DPA1, we will put the tiles in between
.HLACNV.syncTileVar <- function(tile.dt, var.dt, HLACNV.opts) {
    coords <- .HLACNV.getConsts()
    tile.dt[, order.tile := 1:.N] # for QC later
    var.dt [, order.var := 1:.N]   # for QC later
    tile.split <- split(tile.dt, tile.dt$group)
    var.split <- split(var.dt, var.dt$group)

    coords <- coords[coords$gene_name %in% tile.dt$group]

    subdata <- lapply(1:length(coords), function(i) {
        ref.coord      <- coords[i]
        subtile        <- tile.split[[ref.coord$gene_name]]
        subvar         <- var.split[[ref.coord$gene_name]]

        shared.peaks <- intersect(subtile$peak.id, subvar$peak.id)

        if(is.null(shared.peaks)) { ## if no shared peak available
          ## TODO
          subtile$seqnames <- "chr6"
          new.coords       <- floor(seq(start(ref.coord)+HLACNV.opts$tile.shoulder, end(ref.coord)-HLACNV.opts$tile.shoulder, length.out=nrow(subtile)+1))
          subtile$start    <- new.coords[1:nrow(subtile)]+1
          subtile$end      <- new.coords[2:length(new.coords)]
          return(list(syncsubtile = subtile, syncsubvar = subvar))
        }

        subtile <- subtile[peak.id %in% shared.peaks]
        subvar  <- subvar[ peak.id %in% shared.peaks]



        if(is.null(subvar)) { ## if no var available
            ## TODO
            subtile$seqnames <- "chr6"
            new.coords       <- floor(seq(start(ref.coord)+HLACNV.opts$tile.shoulder, end(ref.coord)-HLACNV.opts$tile.shoulder, length.out=nrow(subtile)+1))
            subtile$start    <- new.coords[1:nrow(subtile)]+1
            subtile$end      <- new.coords[2:length(new.coords)]
            return(list(syncsubtile = subtile, syncsubvar = subvar))
        }


        ## distributing tiles in the regions
        subtile$seqnames <- "chr6"
        new.coords       <- floor(seq(start(ref.coord)+HLACNV.opts$tile.shoulder, end(ref.coord)-HLACNV.opts$tile.shoulder, length.out=nrow(subtile)+1))
        subtile$start    <- new.coords[1:nrow(subtile)]+1
        subtile$end      <- new.coords[2:length(new.coords)]

        ##
        subvar.temp <- merge.data.table(subvar[,c("peak.id", "order.var")], subtile[,c("peak.id", "start", "end", "order.tile")], by = "peak.id")
        subvar.temp <- subvar.temp[order(order.tile, order.var)]
        subvar.temp[, ':='(
            new.start = floor(seq(min(start)+HLACNV.opts$var.shoulder, max(end)-HLACNV.opts$var.shoulder, length.out=.N+1))[1:.N]+1,
            new.end   = floor(seq(min(start)+HLACNV.opts$var.shoulder, max(end)-HLACNV.opts$var.shoulder, length.out=.N+1))[2:(.N+1)]
        ), by = peak.id]

        subvar$seqnames  <- "chr6"
        subvar$start     <- (subvar.temp$new.start+subvar.temp$new.end)/2 ## putting in the middle of each tile
        subvar$end       <- (subvar.temp$new.start+subvar.temp$new.end)/2 ## putting in the middle of each tile


        setkey(subtile, seqnames, start, end)
        setkey(subvar, seqnames, start, end)
        new.hits  <- data.table::foverlaps(subtile, subvar, by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"), which = TRUE)

        # if (all(new.hits == hits, na.rm = TRUE)) {print("Sync QC Passed!")}

        return(list(syncsubtile = subtile, syncsubvar = subvar))
    })

    synctile <- rbindlist(lapply(subdata, '[[', 1))
    syncvar  <- rbindlist(lapply(subdata, '[[', 2))

    ## add GRCh38
    synctile <- rbind(tile.split[['GRCh38']], synctile)
    synctile <- synctile[order(order.tile)]

    syncvar <- rbind(var.split[['GRCh38']], syncvar)
    syncvar <- syncvar[order(order.var)]

    syncvar$order.var   <- NULL
    synctile$order.tile <- NULL

    synctile$arm <- factor(synctile$arm, levels = unique(synctile$arm), ordered = TRUE)

    return(list(tile = synctile, var = syncvar))
}



.HLACNV.modelFit <- function(mhla, seg, purity, ploidy, opts) {
    mhla.chr6 <- mhla
    mhla.chr6$tile <- mhla.chr6$tile[seqnames(mhla.chr6$tile) == "chr6"]
    mhla.chr6$var  <- mhla.chr6$var[seqnames(mhla.chr6$var) == "chr6"]
    data <- .HLACNV.opt.data(mhla.chr6, seg, opts)
    fit <- .HLACNV.opt.fit(data, purity, ploidy, opts)
    return(fit)
}



.HLACNV.opt.data <- function(mcnv, seg, opts) {
    var.seg <- findOverlaps(mcnv$var, seg, select="first")
    tile.seg <- findOverlaps(mcnv$tile, seg, select="first")
    if (length(mcnv$var)>0) {
        ## TODO: it is possible that this will return a 0-row table, fix.
        af <- cnvex:::.af.opt.data(mcnv$var, var.seg)
    } else {
        af <- NULL
    }
    lr <- cnvex:::.lr.opt.data(mcnv$tile, tile.seg)
    seg <- cnvex:::.seg.opt.data(seg, mcnv$tile, tile.seg)
    data <- list(lr=lr, af=af, seg=seg, stats=mcnv$stats)
    return(data)
}

.HLACNV.opt.fit <- function(data, p, P, opts) {
    lrC <- cnvex:::.lr.grid.lrC(data$lr, data$stats$sd.lr, opts$opt.max.C)
    lrl.pick <- cnvex:::.llik.lrC.outer(lrC, p, P, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE, opts$max.C)[data$seg]
    afC <- cnvex:::.af.grid.afC(lrl.pick[,.(seg, C)], data$af, p, P, opts$opt.max.C)
    setkey(lrl.pick, seg, C)
    if (!is.null(data$af) && nrow(data$af) > 0) {
        afl.pick <- cnvex:::.llik.afC.outer(afC, opts$opt.p.af.anom, opts$opt.dp.af.max, TRUE)
        setkey(afl.pick, seg, C)
        tmp <- afl.pick[lrl.pick]
        tmp[nC==1, ":="(K=0)]
        fit <- tmp[,.(seg, C=round(C * nC/2), K, lr, tL=i.llik, aL=llik, aD, anom, mse, nlr, naf, len, sC=sC * nC/2)]
    } else {
        fit <- lrl.pick[,.(seg, C=round(C * nC/2), K=NA_integer_, lr, tL=llik, aL=NA_real_, aD, anom=NA_real_,
                           mse=NA_real_, nlr, naf=NA_real_, len, sC=sC * nC/2)]
    }
    return(fit)
}



# work: 1. wgs script
#       2. generating plots for wgs
#       3. get started on Ovarian, pei request
#       4. cnv data for MHC
#       5. Ovarian co-occurance
#       6. Ovarian other plots


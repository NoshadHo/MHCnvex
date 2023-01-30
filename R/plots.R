.HLACNV.getCols <- function() {
    MHCCOLORS <- c(
        GRCh38       = "#b2ccc6",
        `A-HLA`    = "#3fe34b",
        `B-HLA`    = "#777cd2",
        `C-HLA`    = "#bd395c",
        `DPA1-HLA` = "#b8e27d",
        `DPB1-HLA` = "#e855d0",
        `DQA1-HLA` = "#256676",
        `DQB1-HLA` = "#5826a6",
        `DRA-HLA`  = "#429369",
        `DRB1-HLA` = "#7c225f"

    )
    return(MHCCOLORS)
}



.HLACNV.plotZoom <- function(pdata, HLACNV.opts) {
    group.seq <- rle(pdata$group)
    hla.idx   <- cumsum(group.seq$lengths)
    zdata     <- pdata[1:hla.idx[1]]
    for(i in 2:length(hla.idx)) {
        start.idx    <- hla.idx[i-1]+1
        end.idx      <- hla.idx[i]
        region       <- pdata[start.idx:end.idx]
        base.x.start <- zdata[start.idx-1]$x.coord
        if(all(region$group == "GRCh38")) {
            region[, x.coord := base.x.start + cumsum(c(0,diff(x.coord)))]
        } else {
            region[, x.coord := base.x.start + cumsum(c(HLACNV.opts$intplot.zoom.fold, diff(x.coord) * HLACNV.opts$intplot.zoom.fold))]
        }

        zdata <- rbind(zdata, region)
    }
    ggplot(data = zdata)+geom_point(aes(x = x.coord, y = AF, color = group))

    return(zdata)
}


.HLACNV.plotBaf <- function(mhla, gene, p, P, plot.type = c("t.AF", "n.AF"), HLACNV.opts) {
    coords <- .HLACNV.getConsts()
    xlims <- c(min(start(coords))-HLACNV.opts$intplot.left.range*1,
               max(end(coords))  +HLACNV.opts$intplot.right.range*1)
    pdata  <- as.data.table(mhla$var)

    ## set the range
    pdata  <- pdata[start >= min(start(coords))-HLACNV.opts$intplot.left.range &
                        end <= max(end(coords))+HLACNV.opts$intplot.right.range]

    ## set the x.coord
    if (HLACNV.opts$intplot.plot.x.class == "index") {
        pdata[, x.coord := 1:.N]
    } else if (HLACNV.opts$intplot.plot.x.class == "coords") {
        pdata[, x.coord := start]
    }

    ## zooming
    # pdata <- .HLACNV.plotZoom(pdata, HLACNV.opts)

    ## switch coords
    group.seq    <- rle(pdata$group)
    switch.idx   <- cumsum(group.seq$lengths)
    switch.idx = switch.idx+1
    if(HLACNV.opts$intplot.plot.x.class == "coords") {
        switch.coord <- data.table(group = group.seq$values, vline.coord = pdata[switch.idx]$x.coord)
    } else {
        switch.coord <- data.table(group = group.seq$values, vline.coord = pdata[switch.idx]$x.coord)
    }
    switch.coord <- head(switch.coord, -1)
    switch.coord$vline.coord = switch.coord$vline.coord + HLACNV.opts$intplot.zoom.fold/2 ## between two point
    switch.coord[, line.style := "dotted"]
    switch.coord[group == "GRCh38", line.style := "solid"]

    label.coord = data.table(x     = (c(min(pdata$x.coord)-1,switch.coord$vline.coord) + c(switch.coord$vline.coord, max(pdata$x.coord)))/2,
                             y     = Inf,
                             label = c(switch.coord$group, "GRCh38"),
                             # vjust = rep_len(c(4,1),nrow(switch.coord)+1))
                             vjust = 1)
    label.coord[label != "GRCh38", label := tstrsplit(label, "-HLA")[1]]
    label.coord[label == "GRCh38", ':='(label = "chr6", y = -Inf, vjust = -1)]


    grid <- expand.grid(C = unique(gene$C), K = unique(gene$K))
    vars.dip <- .HLACNV.baseLogRatioVars(p, P, c(0,8), diploid = TRUE,Msel = grid$K, Csel = grid$C)
    MCbaf <- vars.dip$MCbaf
    MHCCOLORS <-.HLACNV.getCols()
    # gene.text <- paste(sprintf("%s: C=%s, K=%s",gene$group, gene$C, gene$K),collapse = "\n")
    if (plot.type == "t.AF") {
        pdata$pd <- pdata[,c("t.AF")]

        plt <- ggplot(data=pdata)+
            geom_point(aes(x = x.coord, y = pd, color = group, size = group))+scale_size_manual(values = c(rep(1.5,9), 0.5))+
            geom_hline(yintercept = c(0,1))+
            geom_hline(yintercept = c(0.5), linetype = "dotted", color = "grey")+
            geom_hline(yintercept = MCbaf$baf, linetype = "dotted", color = "grey")+
            # geom_vline(xintercept = switch.coord$vline.coord, color = "black", linetype = switch.coord$line.style)+
            scale_y_continuous(breaks=MCbaf$baf,labels = MCbaf$lab)+
            theme_pubr()+scale_color_manual(values = MHCCOLORS)+
            geom_text_repel(data = label.coord, aes(x = label.coord$x, y = label.coord$y, label = label.coord$label), size = 3,
                            color = "darkred", fontface = "bold")+
            theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")

    } else if (plot.type == "n.AF") {
        pdata$pd <- pdata[,c("n.AF")]

        plt <- ggplot(data=pdata)+
            geom_point(aes(x = x.coord, y = pd, color = group, size = group))+scale_size_manual(values = c(rep(1.5,9), 0.5))+
            geom_hline(yintercept = c(0,1))+
            geom_hline(yintercept = c(0.5), linetype = "dotted", color = "grey")+
            geom_hline(yintercept = MCbaf$baf, linetype = "dotted", color = "grey")+
            # geom_vline(xintercept = switch.coord$vline.coord, color = "black", linetype = switch.coord$line.style)+
            scale_y_continuous(breaks=MCbaf$baf,labels = MCbaf$lab)+
            theme_pubr()+scale_color_manual(values = MHCCOLORS)+
            geom_text_repel(data = label.coord, aes(x = label.coord$x, y = label.coord$y, label = label.coord$label), size = 3,
                            color = "darkred", fontface = "bold")+
            theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")

    }


}


.HLACNV.plotTile <- function(mhla, gene, seg, p, P, HLACNV.opts) {
    coords <- .HLACNV.getConsts()
    xlims <- c(min(start(coords))-HLACNV.opts$intplot.left.range*1,
               max(end(coords))  +HLACNV.opts$intplot.right.range*1)
    pdata  <- as.data.table(mhla$tile)

    ## set the range
    pdata  <- pdata[start >= min(start(coords))-HLACNV.opts$intplot.left.range &
                        end <= max(end(coords))+HLACNV.opts$intplot.right.range]

    ## set the x.coord
    if (HLACNV.opts$intplot.plot.x.class == "index") {
        pdata[, x.coord := 1:.N]
    } else if (HLACNV.opts$intplot.plot.x.class == "coords") {
        pdata[, x.coord := start]
    }

    ## zooming
    # pdata <- .HLACNV.plotZoom(pdata, HLACNV.opts)


    ## switch coords
    group.seq    <- rle(pdata$group)
    switch.idx   <- cumsum(group.seq$lengths)
    switch.idx = switch.idx+1
    if(HLACNV.opts$intplot.plot.x.class == "coords") {
        switch.coord <- data.table(group = group.seq$values, vline.coord = pdata[switch.idx]$x.coord)
    } else {
        switch.coord <- data.table(group = group.seq$values, vline.coord = pdata[switch.idx]$x.coord)
    }
    switch.coord <- head(switch.coord, -1)
    switch.coord$vline.coord = switch.coord$vline.coord + HLACNV.opts$intplot.zoom.fold/2 ## between two point
    switch.coord[, line.style := "dotted"]
    switch.coord[group == "GRCh38", line.style := "solid"]

    label.coord = data.table(x     = (c(min(pdata$x.coord)-1,switch.coord$vline.coord) + c(switch.coord$vline.coord, max(pdata$x.coord)))/2,
                             y     = Inf,
                             label = c(switch.coord$group, "GRCh38"),
                             # vjust = rep_len(c(4,1),nrow(switch.coord)+1))
                             vjust = 1)
    label.coord[label != "GRCh38", label := tstrsplit(label, "-HLA")[1]]
    label.coord[label == "GRCh38", ':='(label = "chr6", y = -Inf, vjust = -1)]


    vars.dip <- .HLACNV.baseLogRatioVars(p, P, c(0,8), diploid = TRUE)
    Clr <- vars.dip$Clr

    pdata$pd <- pdata[,c("lr")]
    MHCCOLORS <-.HLACNV.getCols()

    ## adding segments
    seg.tbl <- as.data.table(seg)
    setkey(pdata, seqnames, start, end)
    hits <- foverlaps(seg.tbl, pdata, by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"), which = TRUE)
    pdata[hits[!is.na(yid)]$yid, seg := hits[!is.na(yid)]$xid]
    plot.seg <- pdata[,.(x.start = min(x.coord), x.end = max(x.coord), y = mean(pd,na.rm=TRUE)), by = seg]

    plt <- ggplot(data=pdata)+
        geom_point(aes(x = x.coord, y = pd, color = group, size = group))+scale_size_manual(values = c(rep(1.5,9), 0.5))+
        geom_hline(yintercept = Clr$lr, color = "grey", linetype = "dotted")+
        geom_segment(data = plot.seg, aes(x = x.start, xend = x.end,y = y, yend = y),colour = 'red')+
        # geom_vline(xintercept = switch.coord$vline.coord, color = "black", linetype = switch.coord$line.style)+
        theme_pubr()+scale_color_manual(values = MHCCOLORS)+
        scale_y_continuous(breaks=Clr$lr,labels = Clr$C)+
        geom_text_repel(data = label.coord, aes(x = label.coord$x, y = label.coord$y, label = label.coord$label), size = 3,
                        color = "darkred", fontface = "bold")+
        theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")



}

.HLACNV.plotLrBaf <- function(mhla, gene, seg, p, P, plot.type = c("t.AF", "n.AF"), HLACNV.opts) {
    mhla.chr6 <- mhla
    mhla.chr6$tile <- mhla.chr6$tile[seqnames(mhla.chr6$tile) == "chr6"]
    mhla.chr6$var  <- mhla.chr6$var[seqnames(mhla.chr6$var) == "chr6"]

    p.baf <- .HLACNV.plotBaf(mhla.chr6, gene, p, P, plot.type = "t.AF", HLACNV.opts)
    p.lr  <- .HLACNV.plotTile(mhla.chr6, gene, seg, p, P, HLACNV.opts)

    ## text annotation
    gene.names <- str_pad(unlist(strsplit(gene$group, "-HLA")), 7, side = "right")
    gene.text  <- paste(sprintf("%s\t: C=%d,\tK=%d,\tHet=%s",gene.names, gene$C, gene$K, gene$hetero),collapse = "\n")
    p.gene     <- ggplot()+geom_label(aes(x = 0, y = 0),label = gene.text, data=gene, hjust = 0)+theme_void()

    (p.lr/p.baf)+inset_element(p.gene, left = 0.6, bottom = 0.6, top = 1.5, right = 1)

}


.HLACNV.baseLogRatioVars <- function(purity, ploidy, range, Msel=0, Csel=1, diploid=TRUE) {
    p <- purity
    D <- (ploidy * p) + 2 * (1-p)
    C <- seq(range[1], range[2], by=1)
    if (diploid == FALSE) {
        Clr <- data.table(C=factor(C), lr=log2((p * C + (1-p) * 1)  / D))
        ymin <- log2((p * min(C) + (1-p) * 1) / D) - 0.25
        ymax <- log2((p * max(C) + (1-p) * 1) / D) + 0.25
    } else {
        Clr <- data.table(C=factor(C), lr=log2((p * C + (1-p) * 2)  / D))
        ymin <- log2((p * min(C) + (1-p) * 2) / D) - 0.25
        ymax <- log2((p * max(C) + (1-p) * 2) / D) + 0.25
    }

    baf <- (p * Msel + 1 * (1-p)) / (p * Csel + 2 * (1-p))
    MCbaf <- data.table(
        baf=c(baf,1-baf),
        lab=rep(paste(Msel,Csel,sep="/"), 2)
    )
    vars <- list(p=p, D=D, C=C, Clr=Clr, ymin=ymin, ymax=ymax, MCbaf=MCbaf)
    return(vars)
}



## TODO: plot for making peaks
# genes <- cov.split[[1]]$hla.type[cov.split[[1]]$hla.type %in% c("A","B","C")]
# a = (cov.mx[, which(cov.split[[1]]$hla.type %in% c("A","B","C"))])
#
# col = rev(RColorBrewer::brewer.pal(11, "Spectral"))
# mean.vals <- colMeans(a, na.rm = TRUE)
# Heatmap(name = "Coverage", a,
#         col = rev(RColorBrewer::brewer.pal(ncol(a), "Spectral")),
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
#         row_title = "Samples", column_title = "Genome Coordinate",
#         raster_quality = 1.5, use_raster = TRUE)
# decorate_heatmap_body("Coverage", {
#     x = (1:ncol(a))/ncol(a)
#     mean.vals.scale <- (mean.vals - min(mean.vals))/(max(mean.vals) - min(mean.vals))
#     grid.lines(x, mean.vals.scale, gp = gpar(lwd = 2, lty = 2))
#     # sapply(peaks.c[peak.cov.PASS==TRUE]$peak.index, function(i) lines(c(x[i]), c(mean.vals.scale[i]), col="Red", lty=2))
#     # sapply(peaks.c[peak.cov.PASS==FALSE]$peak.index, function(i) lines(c(x[i],x[i]), c(y.min, yhat[i]), col="black", lty=2))
# })

# genes   <- pool$gene[pool$gene %in% c("A","B","C")]
# plt.cov <- (pool$coverages[, which(pool$gene %in% c("A","B","C"))])
#
# col = rev(RColorBrewer::brewer.pal(11, "Spectral"))
# mean.vals <- colMeans(plt.cov, na.rm = TRUE)
# Heatmap(name = "Coverage", plt.cov,
#         col = rev(RColorBrewer::brewer.pal(ncol(plt.cov), "Spectral")),
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
#         row_title = "Samples", column_title = "Genome Coordinate",
#         raster_quality = 1.5, use_raster = TRUE)
# decorate_heatmap_body("Coverage", {
#     x = (1:ncol(plt.cov))/ncol(plt.cov)
#     mean.vals.scale <- (mean.vals - min(mean.vals))/(max(mean.vals) - min(mean.vals))
#     grid.lines(x, mean.vals.scale, gp = gpar(lwd = 2, lty = 2))
#     # sapply(peaks.c[peak.cov.PASS==TRUE]$peak.index, function(i) lines(c(x[i]), c(mean.vals.scale[i]), col="Red", lty=2))
#     # sapply(peaks.c[peak.cov.PASS==FALSE]$peak.index, function(i) lines(c(x[i],x[i]), c(y.min, yhat[i]), col="black", lty=2))
# })

#


.HLACNV.plotPool <- function(pool, HLACNV.opts) {

    # genes   <- pool$gene
    # plt.cov <- pool$coverages
    # col = rev(RColorBrewer::brewer.pal(11, "Spectral"))
    # peaks     <- ifelse(is.na(pool$peak.id), FALSE, TRUE)
    # exons     <- ifelse(is.na(pool$ref.exons), FALSE, TRUE)

    genes   <- pool$gene[pool$gene %in% c("A","B","C")]
    plt.cov <- (pool$coverages[, which(pool$gene %in% c("A","B","C"))])
    mean.vals <- colMeans(plt.cov, na.rm = TRUE)
    peaks     <- !is.na(pool$peak.id[which(pool$gene %in% c("A","B","C"))])
    exons     <-  !is.na(pool$ref.exons[which(pool$gene %in% c("A","B","C"))])

    exon.annotation <- HeatmapAnnotation(
        `Coding Region` = exons, `Peak Region` = peaks,

        col = list(`Coding Region` = c("TRUE" = "black", "FALSE" = "white"),
                   `Peak Region`   = c("TRUE" = "darkred", "FALSE" = "white"))
    )

    ht = draw(Heatmap(name = "Coverage", plt.cov,
                 col = rev(RColorBrewer::brewer.pal(ncol(plt.cov), "Spectral")),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                 row_title = "Samples", column_title = "Genome Coordinate",
                 raster_quality = 7, use_raster = TRUE,
                 bottom_annotation = exon.annotation))
    decorate_heatmap_body("Coverage", {
        x = (1:ncol(plt.cov))/ncol(plt.cov)
        mean.vals.scale <- (mean.vals - min(mean.vals))/(max(mean.vals) - min(mean.vals))
        grid.lines(x, mean.vals.scale, gp = gpar(lwd = 2, lty = 2))
        # sapply(peaks.c[peak.cov.PASS==TRUE]$peak.index, function(i) lines(c(x[i]), c(mean.vals.scale[i]), col="Red", lty=2))
        # sapply(peaks.c[peak.cov.PASS==FALSE]$peak.index, function(i) lines(c(x[i],x[i]), c(y.min, yhat[i]), col="black", lty=2))
    })

    # ht <- draw(ht)
    # decorate_heatmap_body("Coverage", {
    #     x = (1:ncol(plt.cov))/ncol(plt.cov)
    #     mean.vals.scale <- (mean.vals - min(mean.vals))/(max(mean.vals) - min(mean.vals))
    #     grid.lines(x, mean.vals.scale, gp = gpar(lwd = 2, lty = 2))
        # sapply(peaks.c[peak.cov.PASS==TRUE]$peak.index, function(i) lines(c(x[i]), c(mean.vals.scale[i]), col="Red", lty=2))
        # sapply(peaks.c[peak.cov.PASS==FALSE]$peak.index, function(i) lines(c(x[i],x[i]), c(y.min, yhat[i]), col="black", lty=2))
    # })

    return(invisible(NULL))
}

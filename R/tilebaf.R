.HLACNV.add.baf <- function(tile.dt, var.dt, opts) {
    ## read AF from var and add BAF to tile
    tile.dt$order <- 1:nrow(tile.dt)
    tile.split <- split(tile.dt, (tile.dt$group))
    var.split  <- split(var.dt, (var.dt$group))

    tile.dt <- lapply(1:length(tile.split), function(i) {
        subtile <- tile.split[[i]]
        subtile[,in.var := FALSE]
        subvar  <- var.split[[subtile$group[1]]]
        if(is.null(subvar)) { ## if no var available
            subtile$baf <- NA_real_
            subtile$baf.depth <- NA_real_
            subtile$baf.n <- NA_real_
            return(subtile)
        }



        if (subtile$group[1] == "GRCh38") {
            subtile <- .HLACNV.mirrorBAF(subtile, subvar, opts)
        } else {
            ## count var on each peak
            count <- subvar[,.(n=.N, start=min(start), end = max(end)),by=peak.id]
            dev.null <- lapply(1:nrow(count), function(j) {
                row <- count[j,]
                subtile[peak.id == row$peak.id, ':='(start = row$start, end = row$end, in.var=TRUE)]
                # return(subtile[peak.id == row$peak.id,])
            })
            ## handling  peaks in tile not in var
            ## preventing duplicated start/end coords
            subtile[in.var == FALSE, ':='(start=start*10000, end=end*10000)]

            subtile <- .HLACNV.mirrorBAF(subtile, subvar, opts)
        }
        return(subtile)
    })
    tile.dt <- rbindlist(tile.dt)
    tile.dt <- tile.dt[order(order)]
    # setkey(tile, order)
    # tile$order <- NULL
    return(tile.dt)
}

.HLACNV.mirrorBAF <- function(tile.tbl, var.tbl, opts){
    tile <- makeGRangesFromDataFrame(tile.tbl, keep.extra.columns = TRUE)
    var  <- makeGRangesFromDataFrame(var.tbl, keep.extra.columns = TRUE)
    hits <- .HLACNV.getHits(var, tile, opts)
    if (length(hits)>0) {
        bad=ifelse(var[queryHits(hits)]$AF < 0.5,
                   round(   var[queryHits(hits)]$AF  * var[queryHits(hits)]$DP),
                   round((1-var[queryHits(hits)]$AF) * var[queryHits(hits)]$DP))
        tmp <- data.table(
            idx=subjectHits(hits),
            bad=bad,
            depth=var[queryHits(hits)]$DP
        )
        setkey(tmp, idx)
        tmp <- tmp[J(1:length(tile))]
        tmp <- tmp[,.(
            baf=sum(bad)/sum(depth),
            depth=sum(depth),
            n=length(na.omit(bad))
        ), by=idx]
        tile.tbl$baf <- cnvex:::.smoothOutliers(tmp$baf, tile, opts)
        tile.tbl$baf.depth <- tmp$depth
        tile.tbl$baf.n <- tmp$n
    } else {
        tile.tbl$baf <- NA_real_
        tile.tbl$baf.depth <- NA_real_
        tile.tbl$baf.n <- NA_real_
    }
    return(tile.tbl)
}

.HLACNV.getHits <- function(snp, tile, opts) {
    if (any(!tile$target)) {
        hits <- findOverlaps(tile, snp)
        hits <- hits[tile[subjectHits(hits)]$target] # prefer assignment to target
        hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
    } else {
        hits <- findOverlaps(snp, tile)
    }
    return(hits)
}

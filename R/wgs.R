.HLACNV.makeWgsCompatible <- function(syncdata, HLACNV.opts) {
  tile          <- syncdata$tile[syncdata$tile$arm == "chr6p"]
  tile.grp.id   <- rleidv(tile$group)
  tile.split    <- split(tile, tile.grp.id)

  hla.groups    <- sort(unique(tile[tile$group != "GRCh38"]$group))
  hla.fixed.slices <- lapply(hla.groups, function(grp) {
    ## find slices
    hla.slice.idx   <- tile.grp.id[which(tile$group == grp)][1]
    left.slice.idx  <- tile.grp.id[max(which(tile.grp.id < hla.slice.idx & tile$group == "GRCh38"))][1]
    right.slice.idx <- tile.grp.id[min(which(tile.grp.id > hla.slice.idx & tile$group == "GRCh38"))][1]

    hla.slice       <- tile.split[[hla.slice.idx]]
    left.slice      <- tile.split[[left.slice.idx]]
    right.slice     <- tile.split[[right.slice.idx]]

    ## find the neighboring regions
    neighborhood    <- .HLACNV.bafFindNeighborhood(hla.slice, HLACNV.opts)

    ## find closest slice to match CN
    sel.slice       <- .HLACNV.bafClosestSlice(hla.slice, left.slice, right.slice)

    ## pick the closest n points
    sel.slice <- sel.slice[start >= min(neighborhood$start) & end <= max(neighborhood$end)]
    if(nrow(sel.slice) == 0) {
        sel.slice <- .HLACNV.bafClosestSlice(hla.slice, left.slice, right.slice)[1]
    }

    hla.tiles <- sel.slice[sample(1:nrow(sel.slice), min(HLACNV.opts$wgs.tile.count, min(nrow(sel.slice), nrow(hla.slice))))]

    ## adjust information
    hla.tiles <- .HLACNV.bafAdjustInfo(hla.tiles, hla.slice, HLACNV.opts)

    return(hla.tiles)
  })
  hla.fixed.slices <- hla.fixed.slices[sapply(hla.fixed.slices, nrow) > 0] ## if no tile is returned, then remove it
  ## replace step
  ## iterationally, replacing tile at place

  for(new.hla.slice in hla.fixed.slices) {
    hla.group <- new.hla.slice$group[1]
    tile.idx  <- which(tile$group == hla.group)
    tile      <- rbind(tile[1:(min(tile.idx)-1)], new.hla.slice, tile[(max(tile.idx)+1):nrow(tile)])
  }

  left.genome   <- syncdata$tile[arm < "chr6p"]
  right.genome  <- syncdata$tile[arm > "chr6p"]
  syncdata$tile <- rbind(left.genome, tile, right.genome)

  return(syncdata)
}

.HLACNV.bafAdjustInfo <- function(hla.tiles, hla.slice, HLACNV.opts) {
  ## replace the coordinations with same number of coordinates from slice
  selected.coords <- hla.slice[sort(sample(1:nrow(hla.slice), nrow(hla.tiles)))]
  hla.tiles[, c("seqnames", "start", "end", "group", "peak.id", "idx", "order")] <-
    selected.coords[, c("seqnames", "start", "end", "group", "peak.id", "idx", "order")]

  return(hla.tiles)
}

.HLACNV.bafFindNeighborhood <- function(hla.slice, HLACNV.opts) {
  left.neighborhood  <- data.table(seqnames = "chr6",
                                   start    = min(hla.slice$start)-HLACNV.opts$wgs.adjacent.range,
                                   end      = min(hla.slice$start)-1 )
  right.neighborhood <- data.table(seqnames = "chr6",
                                   start    = max(hla.slice$start)+1,
                                   end      = max(hla.slice$start)+HLACNV.opts$wgs.adjacent.range )
  neighborhood       <- rbind(left.neighborhood, right.neighborhood)

  return(neighborhood)
}

.HLACNV.bafClosestSlice <- function(hla.slice, left.slice, right.slice) {
  mBAF.hla     <- median(hla.slice$baf,   na.rm=TRUE)
  mBAF.left    <- median(left.slice$baf,  na.rm=TRUE)
  mBAF.right   <- median(right.slice$baf, na.rm=TRUE)

  left.div     <- abs(mBAF.left - mBAF.hla)
  right.div    <- abs(mBAF.right- mBAF.hla)


  ## handle NA cases
  if(is.na(left.div))  {return(right.slice)}
  if(is.na(right.div)) {return(left.slice)}

  ## if everything is good:
  if(left.div > right.div) {return(right.slice)} else {return(left.slice)}
}

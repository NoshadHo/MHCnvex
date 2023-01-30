.HLACNV.segOut <- function(fit, seg, mhla) {
    mcols(seg) <- fit[,.(seg,C,K,sC=round(sC, 4),lr=round(lr, 4),nlr,naf)]
    hits       <- as.data.table(findOverlaps(seg, mhla$tile))
    seg.group  <- as.data.table(mhla$tile)
    seg.group$seg <- NA
    seg.group$seg[hits$subjectHits] <- hits$queryHits
    groups  <- seg.group[,.(group = list(unique(group))),by = seg]
    seg$group  <- groups[!is.na(groups$seg)]$group
    return(seg)
}

.HLACNV.tileOut <- function(fit, seg, mhla) {
    seg       <- .HLACNV.segOut(fit, seg, mhla)
    hits      <- findOverlaps(seg, mhla$tile)
    gene.tile <- mhla$tile[subjectHits(hits)]

    gene.tile$seg <- queryHits(hits)
    gene.tile$C   <- seg[queryHits(hits)]$C
    gene.tile$K   <- seg[queryHits(hits)]$K
    gene.tile$sC  <- seg[queryHits(hits)]$sC

    return(gene.tile)
}

.HLACNV.geneOut <- function(fit, seg, mhla) {
    gene.tile  <- as.data.table(.HLACNV.tileOut(fit,seg,mhla))
    gene.tile  <- gene.tile[group != "GRCh38"]
    g.num      <- gene.tile[,.(n=.N, lr = mean(lr,na.rm=TRUE)),by = .(group,C,K, sC, seg)]
    g.pick     <- g.num[,.SD[which.max(n)],by = group] ## assumption that each gene is no more than 2 segments
    g.pick     <- merge.data.table(g.pick, mhla$haps, by = "group")

    ## find minor, major allele

    return(g.pick)
}

.HLACNV.haploidCnv <- function(basepair, gene, mismatch) {
    haploids <- lapply(1:nrow(gene), function(i) {
        row      <- gene[i,]
        vars     <- mismatch[group == row$group]
        var.summ <- vars[,.(hap1=hap1[1], hap2=hap2[1], BAF.t=mean(BAF.t, na.rm=TRUE)),by = group]
        if (nrow(var.summ) == 0) { ## if gene with no var: pick one randomly
            high.cov.allele <- basepair[hla.type == row$hla.type]

            if(nrow(high.cov.allele) == 0) {
                return(data.table(group=rep(row$group,2), hap = c(row$hap1, row$hap2),
                                  C = c(ifelse(row$hetero==TRUE,row$C-row$K, row$C),
                                        ifelse(row$hetero==TRUE, row$K, NA)),
                                  hetero = row$hetero, call.base="Random"))
            }

            high.cov.allele <- high.cov.allele[,.(t.cov = median(t.cov, na.rm=TRUE)), by = .(seqnames, hla.type)]
            high.cov.allele <- high.cov.allele[order(-t.cov)]
            return(data.table(group=rep(row$group,2), hap = c(high.cov.allele$seqnames[1], high.cov.allele$seqnames[2]),
                              C = c(ifelse(row$hetero==TRUE,row$C-row$K, row$C),
                                    ifelse(row$hetero==TRUE, row$K, NA)),
                              hetero = row$hetero, call.base="Max Cov"))
        }

        if(var.summ$BAF.t > 0.5) { ## baf is hap1/...
            return(data.table(group=rep(var.summ$group,2), hap = c(var.summ$hap1, var.summ$hap2),
                              C = c(ifelse(row$hetero==TRUE,row$C-row$K, row$C),
                                    ifelse(row$hetero==TRUE, row$K, NA)),
                              hetero = row$hetero, call.base="BAF"))
        } else {
            return(data.table(group=rep(var.summ$group,2), hap = c(var.summ$hap2, var.summ$hap1),
                              C = c(ifelse(row$hetero==TRUE,row$C-row$K, row$C),
                                    ifelse(row$hetero==TRUE, row$K, NA)),
                              hetero = row$hetero, call.base="BAF"))
        }
    })
    haploids <- rbindlist(haploids)
    haploids[, major.allele := (C == max(C, na.rm=TRUE)), by= group]
    haploids[C[1] == C[2], major.allele := NA, by= group]
}


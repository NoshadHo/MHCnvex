HLACNV.opts <- list(
    panel.mode="wxs", # "wgs" vs "wxs"
    haplotyper.mode = "polytect",  # polytect vs hapster
    mode="base", ## "exon" or "base"
    # base.lr.mode="peakTile", # if mode="base", then this option activates. Options are "peakTile" and "peakBase"

    tile.width=75,
    tile.shoulder=50,
    var.shoulder=5,
    cores=30,

    ## wgs
    wgs.tile.count = 3,
    wgs.adjacent.range = 100000,

    ## lr
    lr.mode = "pair", ## pair, mean, pool

    ## ref
    ref.gemline.impute.thr = 0.7,

    ## gc
    bias.gc.mode="auto",    ## options: "on", "off" "auto". If "auto", we determine should it get GC correction or not based on lr.var and bias.diff
    bias.gc.adjust=c(0.2,0.8),
    bias.gc.adjust.span=0.5,
    bias.gc.adjust.offset=TRUE,
    bias.gc.tilt.percentile=0.96,

    ## basepair
    basepair.zerofix=TRUE,
    basepair.zerofix.offset=0.01,
    basepair.mapq=0,
    basepair.cov.flags= 1796, # try 266 and 0

    ## exons
    exons.fix.homo=TRUE, ## it will make one allele cov=1 and the other as sum of the two when homozygous (in agreement with assumption of all reads coming from one allele)

    ## baf
    baf.min.cov=80,

    ## opt
    opt.max.C=10,
    opt.p.lr.anom=0.01,

    ## filter
    filter.MHC.blacklist = TRUE,

    filter.exon.max.n.dev=0.15, #max deviation from 0.5
    filter.exon.min.gene.hetero.percent=0.55, # at least 55% of exons should pass the filter to consider them as heterozygous
    filter.exon.min.n.cov.sd=0.75,
    filter.exon.min.read.num=6,

    filter.perBase.max.n.dev=0.15, #max deviation from 0.5
    filter.perBase.min.gene.hetero.percent=0.40, # at least 40% of exons should pass the filter to consider them
    filter.perBase.min.n.cov.sd=-0.5,    # normal should have at least mean(sd)-(-0.5)*sd
    filter.perBase.min.t.cov.sd=-0.05,    # tumor should have at least mean(sd)-(0.5)*sd
    filter.perBase.min.read.num=80,
    filter.perBase.loess.smooth=0.04,
    filter.perBase.peak.window=50,    ## parameter w of peak finding
    filter.perBase.peak.adjacent=75,      ## how many point adjacent to the peak(left+1+right, so 2X+1) should pass the filter
    ## model
    model.max.C=10,
    model.max.sC=10,
    model.p.lr.anom=0.01,
    model.p.af.anom=0.001,
    model.dp.af.max=1e5,


    ## adjustPp
    adjustPp.ploidy.dev=0.4,
    adjustPp.purity.dev=0.25,
    adjustPp.min.purity=0.10,
    adjustPp.min.ploidy=1,
    adjustPp.purity.precision=0.02,
    adjustPp.ploidy.precision=0.02,
    adjustPp.remove.outlier.lr=TRUE,

    ## outlier
    outlier.min.points.n=3,
    outlier.sd.dev=1,

    ## integrated plot
    intplot.left.range=5000000,
    intplot.right.range=5000000,
    intplot.zoom.fold=2, ## distance of HLA coords will be "intplot.zoom.fold" times wider than others
    intplot.normalized.cov=FALSE,
    intplot.plot.x.class = "index", # index:x to be index of data points, coords: x to be the coordinates

    pool.min.cov.normalized = 15
)

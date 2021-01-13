library(Rsamtools)
library(ape)
library(ggplot2)
# Arguments:
# 1) bam file
# 2) GFF3 annotations
args = commandArgs(trailingOnly=TRUE)
gff <- read.gff(args[2])
gff <- gff[, c('seqid', 'type', 'start', 'end')]
# Convert GFF to wide format, as it should be
gff <- reshape(gff, idvar = "seqid", timevar = 'type', direction = "wide")

bamFile <- BamFile(args[1])

parse_gene = function(row) {
  startcodon = data.frame(pos = c(), len = c(), type = c())
  stopcodon = data.frame(pos = c(), len = c(), type = c())
  gene = as.character(row['seqid'])
  end = as.numeric(row['end.UTR5'])
  gr <- GRanges(seqnames = gene, ranges = IRanges(start = 1, end = end))
  params <- ScanBamParam(which = gr, what = scanBamWhat())
  # find reads aligned to gene of interest
  aln <- scanBam(bamFile, param = params)[[paste0(gene, ':1-',end)]]
  
  if(length(aln[['pos']]) > 0) {
    rib <- data.frame(pos=aln[['pos']], len = aln[['qwidth']])
    # filter for reads that span the start codon
    rib <- rib[rib$pos + rib$len > end,]
    rib <- rib[rib$pos >= end - 14 & rib$pos <= end - 10,]
    rib$pos <-  rib$pos - end - 1
    if(length(rib$pos) > 0) {
      rib$type <- 'startcodon'
    }
    startcodon = rib
  }
  
  stop = as.numeric(row['end.CDS'])
  gr <- GRanges(seqnames = gene, ranges = IRanges(start = stop - 50, end = stop))
  params <- ScanBamParam(which = gr, what = scanBamWhat())
  aln <- scanBam(bamFile, param = params)[[paste0(gene, ':', stop - 50, '-', stop)]]
  if(length(aln[['pos']]) > 0) {
    rib <- data.frame(pos=aln[['pos']], len = aln[['qwidth']])
    rib <- rib[rib$pos + rib$len >= stop, ]
    # filter for reads that span the start codon
    rib$pos = rib$pos + rib$len - stop
    rib <- rib[rib$pos >= 9 & rib$pos <= 13,]
    if(length(rib$pos) > 0) {
      rib$type <- 'stopcodon'
    }
    stopcodon = rib
  }
  
  return(rbind(startcodon, stopcodon))
}

rib = do.call("rbind", apply(gff, 1, FUN = parse_gene))
start = rib[rib$type == 'startcodon',c('pos','len')]
stop = rib[rib$type == 'stopcodon',c('pos','len')]
start = as.data.frame(table(start))
ggplot(start, aes(x=pos, y=len)) + geom_tile(aes(fill = Freq)) + scale_fill_gradient(low="white", high="blue") + theme_minimal()
stop = as.data.frame(table(stop))
ggplot(stop, aes(x=pos, y=len)) + geom_tile(aes(fill = Freq)) + scale_fill_gradient(low="white", high="blue") + theme_minimal()
# aln <- scanBam(bamFile)[[1]]

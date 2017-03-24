#### Determine exonic gene lengths to calculate FPKM.

setwd("~/Shuo_Katja")

# Import the GTF-file that you have also used as input for htseq-count
library(GenomicFeatures)
library(rtracklayer) # version 1.25.1 needed!!


gtf <- "/srv/gsfs0/projects/brunet/Shared/Genomes/C.elegans/ftp.ensembl.org/pub/release-84/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.84.gtf"

txdb <- makeTxDbFromGFF(file = gtf,
                        format="gtf")
gtf <- read.table(gtf, sep = "\t")
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gtf$gene <- gsub("gene_id ", "", sapply(strsplit(as.character(gtf$attribute), ";"), function(x) x[1]))

## Collect the exons per gene id
exons <- gtf[gtf$feature == "exon",]
exons.gr <- GRanges(seqnames=exons$seqname,
                    ranges = IRanges(start=exons$start, end=exons$end),
                    strand=exons$strand,
                    gene = as.factor(exons$gene),
                    frame = exons$frame,
                    attribute = exons$attribute
                    )

## GRangesList: A GRanges object with all exons per gene.
exons.gr.per.gene <- split(exons.gr, f = exons.gr$gene)

## for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
gene.sizes <- sapply(exons.gr.per.gene,
                            function(x){
                                sum(width(reduce(x)))
                            }
                            )
save(gene.sizes,
     file = "Data/gene_sizes.RData")

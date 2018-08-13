#!/usr/bin/env Rscript
#Check this GOOOOOD EXPLANATION OF REPLICATES IN RNASEQ FIRST http://rnaseq.uoregon.edu/
#DeSeq in R
setwd("D:/Google Drive/PhD Research proposal/Objective 4-ToxinDet/results/RNAseqMap")
sampleTable <- read.csv("sample_table.csv",row.names=1)
filenames<-(paste(sampleTable$Run)) #Calling filenames from the sample file $run column
# filenames<-as.character(c("aligned-AF70-CORN.bam", "aligned-AF70-PDB.bam")) #Another option to call for the filenames
file.exists(filenames)

library("Rsamtools")
# bamfiles <- BamFileList(filenames, yieldSize=2000000) #Loading BAM files produced by samtools from the STAR mapping
bamfiles <- BamFileList(filenames) #Loading BAM files produced by samtools from the STAR mapping

seqinfo(bamfiles[1]) #Checking chromosome or contig names


library("GenomicFeatures")
txdb <- makeTxDbFromGFF("af70_20130423.gff3", format="gff3", circ_seqs=character())
ebg <- exonsBy(txdb, by="gene")

#Now we are going to count
library("GenomicAlignments")
library("BiocParallel")
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE,
                        fragments=FALSE)
						
head(assay(se), 3)
colData(se)
colData(se) <- DataFrame(sampleTable)
se$dex
se$dex <- relevel(se$dex, "untrt")
#Quickly check the millions of fragments that uniquely aligned to the genes 
round( colSums(assay(se)) / 1e6, 1 )

#Now using DeSeq
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ dex)

countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)

#construct the DESeqDataSet object from the matrix of counts and the sample information table, we use:

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ dex)
								  
nrow(dds)


plotSparsity(dds, normalized = TRUE)
#rld <- rlog(dds, blind=TRUE) #Did not work and suggested using varianceStabilizingTransformation
rld <- varianceStabilizingTransformation(dds, blind=TRUE)

head(assay(rld), 3)

#Plotting log2 vs the varianceStabilizingTransformation
#Shown are scatterplots using the log2 transform of normalized
#counts (left side) and using the rlog (right side).
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)
	 
#Sample distances
sampleDists <- dist( t( assay(rld) ) )
sampleDists
#Heatmap of sample-to-sample distances using the rlog-transformed values.
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Heatmap of sample-to-sample distances using the Poisson Distance.
	 
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)

#PCA plot Usually used when doing different samples, in this case I am using only two samples and it is not worth it
plotPCA(rld, intgroup = "dex")


#Differential expression analysis

dds <- DESeq(dds)
res <- results(dds)


mcols(res, use.names=TRUE)
summary(res) ##HERE IT WILL SHOW HOW MANY GENES WERE UPREGULATED DUE TO THE PRESENCE OF CORN AND DOWNREGULATED, I WILL USE THESE GENES TO GENERATE E-PROBES

#bEING MORE STRICT AND LOWERING THE FALSE DISCOVERY RATE THRESHOLD 
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

#If we want to raise the log2 fold change threshold, 
#so that we test for genes that show more substantial changes due to treatment
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)


####NOW COMES THE BEST PLOTTING RESULTS

#TOP GENE 
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))

#An MA-plot of changes induced by treatment

plotMA(res, ylim=c(-7,7))
png(filename="MAPlot_GeneChanges.png",bg = "white", width=6000, height=4000, res=500)
plotMA(res, ylim=c(-7,7))
graphics.off()


plotMA(resLFC1, ylim=c(-5,5))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


#Top genes : sort by pvalue : from here: http://www.sthda.com/english/wiki/rna-sequencing-data-analysis-counting-normalization-and-differential-expression
resSort <- res[order(res$pvalue),]
head(resSort)
write.csv(resSort, file="TopGenes-orderedbyPvalue.csv")
resSortlog2 <- res[order(res$log2FoldChange),]
write.csv(resSortlog2, file="TopGenes-orderedbylog2Fold.csv")

#!/usr/bin/env Rscript

#=========================================================================================================================================================
#MAPPING ILLUMINA RNA-SEQ READS TO THE S MINOR GENOME AND DATA ANALYSIS IN R
#Source: https://www.bioconductor.org/help/course-materials/2014/SeattleOct2014/B02.1.1_RNASeqLab.html
#=========================================================================================================================================================

setwd("D:/Google Drive/PhD Research proposal/Objective 4-ToxinDet/results/RNAseqMap")

#------------------------------------------------------------------------------------------------------------------------
#PART 1:  Getting the information for the example
#------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------
#1.1 Sample table
#-------------------------------------------------
sampleTable <- read.csv("sample_table.csv",row.names=1)
filenames<-(paste(sampleTable$Run)) #Calling filenames from the sample file $run column
file.exists(filenames) #Checking if file exists

#-------------------------------------------------
#1.2 Loading BAM files
#-------------------------------------------------
library("Rsamtools")
bamfiles <- BamFileList(filenames)
seqinfo(bamfiles[1])

#-------------------------------------------------
#1.3 Read in the gene model which will be used for counting reads
#-------------------------------------------------
library("GenomicFeatures")
txdb <- makeTxDbFromGFF("af70_20130423.gff3", format="gff3")
(genes <- exonsBy(txdb, by="gene"))

#-------------------------------------------------
#1.3 Counting mapped reads
#-------------------------------------------------
library("GenomicAlignments")
library("BiocParallel")
se <- summarizeOverlaps(features=genes, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE,
                        fragments=FALSE)
se
head(assay(se))
colSums(assay(se))
colData(se)
rowRanges(se)
str(metadata(rowRanges(se)))
(colData(se) <- DataFrame(sampleTable))

#-------------------------------------------------
#2. Running DeSeq2
#-------------------------------------------------
countdata <- assay(se)
head(countdata)
coldata <- colData(se)

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ dex)
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ dex))
								 
								 
# The function rlog returns a SummarizedExperiment object which contains the rlog-transformed values in its assay slot:

rld <- rlog(dds)
head(assay(rld))

# To show the effect of the transformation, we plot the first sample against the second, 
# first simply using the log2 function (after adding 1, to avoid taking the log of zero), 
# and then using the rlog-transformed values. For the log2 method, we need estimate size 
# factors to account for sequencing depth (this is done automatically for the rlog method).

par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )

#--------------------------------------------------------------------------------------------
#2.1 Differential expression analysis	
#--------------------------------------------------------------------------------------------

#Making sure that  untrt is the first level in the dex factor
dds$dex <- relevel(dds$dex, "untrt")
dds <- DESeq(dds) #Running DESeq

#--------------------------------------------------------------------------------------------
# 2.3 Building the results table	
#--------------------------------------------------------------------------------------------
(res <- results(dds))
mcols(res, use.names=TRUE)
summary(res)

#--------------------------------------------------------------------------------------------
# 3. MULTIPLE TESTINGS	
#--------------------------------------------------------------------------------------------
#How many genes genes with a p value below 0.05 among 11534 genes, for which the test succeeded in reporting a p value
sum(res$pvalue < 0.05, na.rm=TRUE)
#Total Genes tested
sum(!is.na(res$pvalue))
#all genes with an adjusted p value below 10 as significant. How many such genes are there?
sum(res$padj < 0.1, na.rm=TRUE)

#................................................................
# We subset the results table to these genes and then sort it
# by the log2 fold change estimate to get the significant genes
# with the strongest down-regulation.
#................................................................
resSig <- subset(res, padj < 0.1)
#Genes with the strongest DOWN-regulation
head(resSig[ order( resSig$log2FoldChange ), ])
write.csv(resSig, file="Genes-strongestDOWNregulation.csv")
#Genes with the strongest UP-regulation
head(resSig[ order( -resSig$log2FoldChange ), ])
UPReg <-resSig[ order( -resSig$log2FoldChange ), ]
write.csv(UPReg, file="Genes-strongestUPregulation.csv")

#--------------------------------------------------------------------------------------------
# 4. DIAGNOSTIC PLOTS	
#--------------------------------------------------------------------------------------------
#......................................................................
# A quick way to visualize the counts for a particular gene is to use
# the plotCounts function, which takes as arguments the DESeqDataSet,
# a gene name, and the group over which to plot the counts.
#......................................................................
 
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))

#...........................................................................
# An “MA-plot” provides a useful overview for an experiment with a two-group
 # comparison. The log2 fold change for a particular comparison is plotted on
 # the y-axis and the average of the counts normalized by size factor is shown
 # on the x-axis (“M” for minus, because a log ratio is equal to log minus log,
 # and “A” for average).
#...............................................................................
plotMA(res, ylim=c(-8,8),colNonSig = "gray", alpha = 0.1)

#Plotting with alpha = 0.1 admits a false discovery rate of 10%, I can add to 20% and will get more red dots
svg(filename = "A.flavus-MAplot.svg")
plotMA(res, ylim=c(-8,8),alpha = 0.1, cex=0.7)
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()


#The function plotDispEsts visualizes DESeq2's dispersion estimates:	 
plotDispEsts(dds)	 

#Another useful diagnostic plot is the histogram of the p values.
hist(res$pvalue, breaks=20, col="grey50", border="white")
#This plot becomes a bit smoother by excluding genes with very small counts:
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")


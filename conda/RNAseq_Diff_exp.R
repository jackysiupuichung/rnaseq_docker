library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

seqdata <- read.delim("GSE60450_LactationGenewiseCounts.txt")
# Read the sample information into R
sampleinfo <- read.delim("SampleInfo.txt", stringsAsFactors = TRUE)

# Remove first two columns from seqdata
countdata <- seqdata[, -(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[, 1]

# shorten colnames to contain only the relevant information about each sample.
colnames(countdata) <- substr(colnames(countdata), start=1, stop=7)

table(colnames(countdata)==sampleinfo$SampleName)

#create a DGEList object, an object used by edgeR to store count data
y <- DGEList(countdata)

#Add the group information into the DGEList
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group <- factor(group)
y$samples$group <- group

##Adding annotation
ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
#confirm the ENTREZID column matches exactly to our y$counts rownames
table(ann$ENTREZID==rownames(y$counts))

##filter out lowly expressed genes
#use cpm function to generate the CPM values and then filter. this normalise for the different sequencing depths for each sample
myCPM <- cpm(countdata)
thresh <- myCPM > 0.5
#keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
#visualise effect of threshold
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5)

##quality control
logcounts <- cpm(y,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#MDSplot, visualisation of a principle components analysis
#ideal scenario: the greatest sources of variation in the data are the treatments/groups we are interested in. It is also an incredibly useful tool for quality control and checking for outliers.
sampleinfo <- read.delim("SampleInfo_Corrected.txt", stringsAsFactors = TRUE)
group <- factor(paste(sampleinfo$CellType,sampleinfo$Status,sep="."))
y$samples$group <- group
par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","black")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

#interactive MDSplot using Glimma package
labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
glMDSPlot(y, labels=labels, groups=group, folder="mds")


# Install packages
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("clusterProfile")

# Load libraries
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ChIPseeker)

# Set working directory
setwd("/home/nuria/Documents/PGB/Chip/Olig2")

############################## Peaks Distribution ############################## 

# Chip profiling
peak <- readPeakFile("Olig2_summits.bed")
covplot(peak, weightCol="V5")

#################### Binding Preference in Genomic Elements #################### 

# Load data
samplefiles <- list(Olig2 ="Peaks.resized.bed")

# Assign annotation database: Mouse genome mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Get annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 1000), verbose=FALSE)

#Plots
plotAnnoPie(peakAnnoList[["Olig2"]])
plotDistToTSS(peakAnnoList, title="DISTRIBUTION OF TRANSCRIPTION FACTOR-BINDING
LOCI RELATIVE TO TSS", face ="bold")

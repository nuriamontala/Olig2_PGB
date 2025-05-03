# Load libraries
library(dplyr)
library(ggplot2)
library(ReactomePA)

setwd("~/Documents/PGB/Chip/Olig2/Chipseeker")
samplefiles <- list(Olig2 = "Peaks.resized.bed")
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 1000), verbose=FALSE)

# Extract annotation data for Olig2 peaks from the peakAnnoList object
Olig2_annot <- as.data.frame(peakAnnoList[["Olig2"]]@anno)

# Extract unique gene IDs (ENTREZ IDs) from the Olig2 annotation data
entrezids <- unique(Olig2_annot$geneId)

# Map ENTREZ IDs to their corresponding gene symbols using the org.Mm.eg.db
database
entrez2gene <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = entrezids,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENTREZID"
)

# Reorder the mapping to match the order of the ENTREZ IDs in the input
entrez2gene <- entrez2gene[match(entrezids, entrez2gene$ENTREZID), c("ENTREZID", "SYMBOL")]
# Rename columns for clarity: "ENTREZID" -> "entrez", "SYMBOL" -> "symbol"
colnames(entrez2gene) <- c("entrez", "symbol")

#Match the gene IDs in Olig2_annot with the mapped ENTREZ-to-symbol data
m <- match(Olig2_annot$geneId, entrez2gene$entrez)

# Add the corresponding gene symbols as a new column to the Olig2_annot data
Olig2_annot <- cbind(Olig2_annot[, 1:14], geneSymbol = entrez2gene$symbol[m])

# Save the updated Olig2 annotation data (with gene symbols) to a text file
write.table(Olig2_annot, file = "Olig2_annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Perform Gene Ontology (GO) enrichment analysis on the ENTREZ IDs
ego <- enrichGO(
  gene = entrezids,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Convert the GO enrichment results to a data frame for easier manipulation
cluster_summary <- as.data.frame(ego)

# Save the enrichment results to a CSV file
write.csv(cluster_summary, "clusterProfiler_Olig2.csv")

# Round the p-value and q-value in the enrichment results for better readability
ego@result <- ego@result %>%
  mutate(pvalue = round(pvalue, 4), qvalue = round(qvalue, 4))

# Create a dot plot to visualize the GO enrichment results
dotplot(ego, showCategory = 25) +
  ggtitle("FUNCTIONAL ENRICHMENT ANALYSIS") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
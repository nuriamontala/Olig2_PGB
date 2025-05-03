# Load libraries
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(pheatmap)

# Set working directory
setwd("~/Documents/PGB/Chip/Olig2")

# Load OLIG2 binding sites
olig2_peaks <- read.table("Peaks.resized.bed", header = FALSE)
olig2_gr <- GRanges(seqnames = olig2_peaks$V1,
                    ranges = IRanges(start = olig2_peaks$V2, end = olig2_peaks$V3))
# Load other TFs
nkx2_2_peaks <- read.table("5_NKX2-2.bed", header = FALSE)
gli3_peaks <- read.table("6_GLI3.bed", header = FALSE)
NKX6_1_peaks <- read.table("7_NKX6-1.bed", header = FALSE)

# Convert other TFs to GRanges
nkx2_2_gr <- GRanges(seqnames = nkx2_2_peaks$V1, ranges = IRanges(start = nkx2_2_peaks$V2, end = nkx2_2_peaks$V3))
gli3_gr <- GRanges(seqnames = gli3_peaks$V1, ranges = IRanges(start = gli3_peaks$V2, end = gli3_peaks$V3))
nkx6_1_gr <- GRanges(seqnames = NKX6_1_peaks$V1, ranges = IRanges(start = NKX6_1_peaks$V2, end = NKX6_1_peaks$V3))

# Calculate Jaccard indices
jaccard_results <- data.frame(
  TF = c("NKX2-2", "GLI3","NKX6_1"),
  Jaccard_Index = c(
    calculate_jaccard(olig2_gr, nkx2_2_gr),
    calculate_jaccard(olig2_gr, gli3_gr),
    calculate_jaccard(olig2_gr, nkx6_1_gr)
  )
)

study_tfs <- list(OLIG2 = olig2_gr, NKX2_2 = nkx2_2_gr, GLI3 = gli3_gr, NKX6_1 = nkx6_1_gr)

# Initialize matrix to store Jaccard indices
jaccard_matrix <- matrix(nrow = length(study_tfs), ncol = length(study_tfs),
                         dimnames = list(names(study_tfs), names(study_tfs)))

# Fill in Jaccard index values for each pair
jaccard_matrix <- matrix(nrow = length(study_tfs), ncol = length(study_tfs),
                         dimnames = list(names(study_tfs), names(study_tfs)))

# Calculate pairwise Jaccard Index for the heatmap
for (i in seq_along(study_tfs)) {
  for (j in seq_along(study_tfs)) {
    jaccard_matrix[i, j] <- calculate_jaccard(study_tfs[[i]], study_tfs[[j]])
  }
}

# Plot heatmap
pheatmap(
  as.matrix(jaccard_matrix),
  color = colorRampPalette(c("white", "#b399c8"))(1000), # Enhanced color palette
  main = "Pairwise Jaccard Index Between OLIG2 and Other TFs",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 12, # Increased font size for better readability
  fontsize_row = 10, # Larger font for row labels
  fontsize_col = 10, # Larger font for column labels
  legend = TRUE, # Show legend
  legend_breaks = seq(0, 1, by = 0.2), # Adjust legend scale
  legend_labels = seq(0, 1, by = 0.2),
  border_color = "grey90", # Subtle borders for cells
  angle_col = 45, # Rotate column labels for readability
  treeheight_row = 30, # Adjust height of row dendrogram
  treeheight_col = 30, # Adjust height of column dendrogram
  display_numbers = TRUE, # Show values in each cell
  number_format = "%.2f", # Format numbers to 2 decimal places
  number_color = "black", # Color of the numbers displayed
  annotation_legend = TRUE # Show annotation legend if used
)
##################### Dot matrices ################################
# Load libraries and set working directory.
library(seqinr)
library(magrittr)
setwd("/home/marti/MBHS/PGB/Project/reviewed_tree/fastas_olig")

# Read all Olig sequences that will be plotted. Olig1/2/3 in human and other Olig2 orthologs.
hs_olig1 <- read.fasta(file = "H._sapiens_Olig1.fasta",seqtype = "AA")[[1]] 
hs_olig2 <- read.fasta(file = "H._sapiens_Olig2.fasta",seqtype = "AA")[[1]] 
hs_olig3 <- read.fasta(file = "H._sapiens_Olig3.fasta",seqtype = "AA")[[1]]

mm_olig2 <- read.fasta(file = "M._musculus_Olig2.fasta",seqtype = "AA")[[1]]
dr_olig2 <- read.fasta(file = "D._rerio_Olig2.fasta",seqtype = "AA")[[1]]
xt_olig2 <- read.fasta(file = "X._tropicalis_Olig2.fasta",seqtype = "AA")[[1]]

# Make a list of the groups of sequences to iterate with.
hs_oligs <- list(hs_olig1, hs_olig2, hs_olig3)

species_olig2 <- list(
  hs_olig2,  # Homo sapiens OLIG2
  mm_olig2,  # Mus musculus OLIG2
  xt_olig2,   # Xenopus tropicalis OLIG2
  dr_olig2  # Danio rerio OLIG2
)
# Try the dotPlot function
#par(mfrow=c(1, 1))
#dotPlot(hs_olig1, ac_olig2, wsize = 15, nmatch = 10)

# Create a vector with the species name for axis labels. 
species_names <- c("H. sapiens", "M. musculus", "D. rerio", "X. tropicalis")
# Create a vector with the names of human Oligs
olig_names <- c("OLIG1 H. sapiens", "OLIG2 H. sapiens", "OLIG3 H. sapiens")


# Configure the grid to plot
par(mfrow = c(4, 3), mar = c(2, 2, 2, 2), oma = c(5, 5, 5, 5))


# Generate the grid of dotPlots with a for loop

for (i in 1:length(species_olig2)) {  # Rows: orthologs
  for (j in 1:length(hs_oligs)) {     # Columns: paralogs in human
    # Create a dotPlot for each combination of sequences from the two lists
    dotPlot(hs_oligs[[j]], species_olig2[[i]], wsize = 15, nmatch = 10)
    
    # Add labels at first column (left of the grid)
    if (j == 1) {
      mtext(species_names[i], side = 2, line = 1, las = 2, cex = 0.8, col = "brown")
    }
    
    # Add labels at the first row (top of the grid)
    if (i == 1) {
      mtext(olig_names[j], side = 3, line = 2, las = 1, cex = 1.2, col = "darkgreen")
    }
    
    # Add a grid inside each dotPlot
    abline(h = seq(0, length(species_olig2[[i]]), by = 50), col = "lightgray")
    abline(v = seq(0, length(hs_oligs[[j]]), by = 50), col = "lightgray")
  }
}

##################### Percent identyty matrix ################################
# Load library for plotting heatmaps
library(pheatmap)

# Read the identity matrix provided by MUSCLE and set as.numeric()
pim <- read.table("/home/marti/MBHS/PGB/Project/reviewed_tree/identity_matrix_corrected.pim", header = FALSE, row.names = 1)
pim[] <- sapply(pim, as.numeric)
colnames(pim) <- rownames(pim)
pim_matrix <- as.matrix(pim)

# Plot the identity matrix as a heatmap in blue gradient
pheatmap(pim_matrix, 
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Percent Identity Matrix", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE
         )



setwd("/home/alexiosgiannoulas/PGB/Project/scRNAseq2")
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(readr)
library(ggplot2)
#Following code can be used with more than one tissue type provided by Tabula Muris, by 
#editing the "grep" step
FACS_files = list.files("FACS/", full.names = TRUE)
FACS_files <- grep("Brain" ,FACS_files, value = TRUE)
raw.data.list <- list()
for (file in FACS_files){
  raw.data <- read.csv(file, row.names = 1)
  raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE) #Very important sparse=TRUE 
  raw.data.list <- append(raw.data.list, raw.data)
}
raw.data <- do.call(cbind, raw.data.list)
#The matrix does not have 0 but "." to save disc space. called "Sparse matrix"

meta.data <- read.csv("metadata_FACS.csv")
#We need to convert this into a metadata table. In the metadata the plate.barcode is the ID.
#We append the names of the plates to the names of the rawdatacols. 
#The colnames of raw.data are the names of the cells.
plates <- str_split(colnames(raw.data),"[.]", simplify = TRUE)[,2]
rownames(meta.data) <- meta.data$plate.barcode
cell.meta.data <- meta.data[plates,]
rownames(cell.meta.data) <- colnames(raw.data) #IMPORTANT the rownames of metadata have to be the same as the colnames of raw.data
#Sequencing preperation included the use of spike-ins
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)#these are spike-ins
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
mt <- grep(pattern = "^Mt", x = rownames(x = raw.data), value = FALSE)
#Because spike ins were not used in the analysis, they are removed from the dataset.
raw.data <- raw.data[-ercc.index,]
#CREATE seurat object
tiss <- CreateSeuratObject(counts = raw.data)
tiss <- AddMetaData(object = tiss, cell.meta.data)
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
#Count percentage of counts originating from a set of features 
tiss[["percent.mt"]] <- PercentageFeatureSet(tiss, pattern = "^Mt")
BeforeQC <- VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA", "percent.ercc", "percent.mt"), ncol = 4) +
  theme(
    plot.title = element_text(size = 14)
  )
ggsave("1beforeQC.png", plot = BeforeQC, width = 8, height = 5, dpi = 300)
tiss <- subset(tiss, subset = nFeature_RNA >500 & nCount_RNA > 50000 & percent.mt < 5)
#From 10561 cells we went to 7460

AfterQC <- VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  theme(
    plot.title = element_text(size = 14)
  )
ggsave("2afterQC.png", plot = AfterQC, width = 6, height = 4, dpi = 300)

tiss <- NormalizeData(tiss, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
tiss <- FindVariableFeatures(object = tiss)
top10 <- head(VariableFeatures(tiss), 10)
plot1 <- VariableFeaturePlot(tiss)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
ggsave("3VariableFeatures.png", plot = plot2, width = 10, height = 7, dpi = 700, bg="white")

#Scaling
tiss <- ScaleData(object = tiss)
#PCA 
tiss <- RunPCA(tiss, features = VariableFeatures(object = tiss))
#Which genes contribute to the variance explained by PC_1 and PC_2
PC1_2_Genes <- VizDimLoadings(tiss, dims = 1:2, reduction = "pca") #Olig1 on PC1
ggsave("4PC1PC2Genes.png", plot = PC1_2_Genes, width = 10, height = 7, dpi = 700, bg="white")

DimPlot(tiss, reduction = "pca") + NoLegend()
Elbow <- ElbowPlot(tiss, ndims = 50)
ggsave("5elbow.png", plot = Elbow , width = 10, height = 7, dpi = 700, bg="white")

#Clustering
tiss <- FindNeighbors(tiss, dims = 1:40)
tiss <- FindClusters(tiss, resolution = 0.5)
#UMAP
tiss <- RunUMAP(tiss, dims = 1:40)
#Plot UMAP for clusters found in previous step
umapclusters <- DimPlot(tiss, reduction = "umap", label = TRUE) 
ggsave("6FirstUMAPclust.png", plot = umapclusters, width = 6, height = 4, dpi = 300)
#Plot UMAP separating by tissue, as separated by FACS sorting
umaptiss <- DimPlot(tiss, reduction = "umap", group.by = 'tissue') #Clear separation between Neurons and Microglia
ggsave("7FirstUMAP.png", plot = umaptiss, width = 6, height = 4, dpi = 300)
#Plot UMAP separating by subtissue
DimPlot(tiss, reduction = "umap", group.by = 'subtissue') #Not clear separation between subtissues

#Annotation using singleR
library(celldex)
library(SingleR)
library(RColorBrewer)
library(SingleCellExperiment)
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
#Transformation of seurat object into sce object, to be compatible with SingleR
sce <- as.SingleCellExperiment(DietSeurat(tiss))
brainmain <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)
brainfine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)
tiss@meta.data$brainmain <- brainmain$pruned.labels
tiss@meta.data$brainfine <- brainfine$pruned.labels
# Define non-brain cell types 
non_brain_cell_types <- c("Adipocytes", "B cells", "Cardiomyocytes", "Erythrocytes", 
                          "Fibroblasts", "Fibroblasts activated", "Fibroblasts senescent", 
                          "Granulocytes", "Hepatocytes", "Monocytes", "NK cells", "T cells")

# Exclude non-brain cell types to avoid false annotation, due to similar transcriptomic profiles
tiss <- tiss[, !(tiss@meta.data$brainfine %in% non_brain_cell_types)]
#Plot UMAP labeled by cell type  
tiss <- SetIdent(tiss, value = "brainfine")
UMAPcelltypes <- DimPlot(tiss, label = T , repel = T, label.size = 4) 
ggsave("8UMAPcelltype.png", plot = UMAPcelltypes, width = 10, height = 7, dpi = 700)
#Plot UMAP showing expression of Olig2
Olig2UMAP <- FeaturePlot(tiss, features = "Olig2")
ggsave("9Olig2UMAP.png", plot = Olig2UMAP, width = 10, height = 7, dpi = 700)
#Save the above 2 UMAPS to compare side by side
ggsave("10Olig2UMAPvsGeneral.png", plot = UMAPcelltypes + Olig2UMAP  , width = 13, height = 8, dpi = 700)

#Plot expression of Olig2 across cell types
features <- c("Olig2")
color_palette <- c("turquoise", "pink", "magenta4", "magenta", "seagreen", "purple", 
                   "violetred1", "cyan", "aquamarine", "deeppink4", "violet", "turquoise4","hotpink")
Olig2Expression <- VlnPlot(tiss, features = features,cols =color_palette )
ggsave("11Olig2violin.png", plot = Olig2Expression, width = 10, height = 7, dpi = 700)

#Save which clusters have each cell type
table_cluster_cellt <-  table(tiss$seurat_clusters, tiss$brainfine)
write.csv(table_cluster_cellt, "clusters_celltypes.csv", row.names = TRUE)
#Find marker genes for cluster 1 which corresponds to oligodendrocytes
tiss <- SetIdent(tiss, value = "seurat_clusters")
tiss.markers.olig<- FindMarkers(tiss, ident.1 = 1, min.pct = 0.25)
write.csv(tiss.markers.olig, "MarkersOligodendrocytes.csv", row.names = T)
#Plot expression of Olig2, Olig1, 5 oligodendeocyte markers, and 1 Microglia marker vs cell types
tiss@meta.data <- tiss@meta.data %>%
  filter(!is.na(brainfine))
#Finding Microglia markers
tiss <- SetIdent(tiss, value = "brainfine")
tiss.markers.Microglia<- FindMarkers(tiss, ident.1 = "Microglia", min.pct = 0.25)
write.csv(tiss.markers.Microglia, "MarkersOligodendrocytes.csv", row.names = T)
#The size of the dot indicates the proportion of cells within the cluster expressing the gene, 
#while the color intensity represents the average expression level of the gene within those cells.
features1<- c("Olig2","Olig1","Ppp1r14a","Ermn","Tmem88b","Aspa","Gpr37","Siglech")
#Siglech known Microglia marker, also found with "FindMarkers", Aspa known Oligodendrocyte marker
colors <- c("turquoise","magenta")
dotplot_brainfine <- DotPlot(tiss, features = features1, group.by = "cell_ontology_class", cols = colors) + RotatedAxis()
ggsave("12DotPlot_brainfine.png", plot = dotplot_brainfine, width = 10, height = 7, ,bg="white", dpi = 700)

#### Annotation according to Tabula Muris
anno <- read_csv("annotations_FACS.csv")
# Assign cell names to the metadata
tiss@meta.data$cell <- rownames(tiss@meta.data)
# Join the metadata with the annotation
meta2 <- tiss@meta.data %>% left_join(anno[, c(1, 3)], by = 'cell')
# Add the cell ontology class metadata
tiss <- AddMetaData(object = tiss, metadata = meta2$cell_ontology_class, col.name = "cell_ontology_class")
# Replace NA values with "unknown"
tiss@meta.data$cell_ontology_class[is.na(tiss@meta.data$cell_ontology_class)] <- "unknown"
# Convert to factor
tiss$cell_ontology_class <- as.factor(tiss$cell_ontology_class)

#Plot UMAP with Tabula Muris annotations
UMAPcelltype2 <- DimPlot(tiss, reduction = "umap", group.by = 'cell_ontology_class')
#Save it together with UMAP-expressionOlig2
ggsave("13Olig2UMAPvsGeneral2.png", plot = UMAPcelltype2 + Olig2UMAP  , width = 13, height = 8, dpi = 700)
#Save it together with UMAP-celltype acording to SingleR
UMAPcelltype1 <- DimPlot(tiss, reduction = "umap", group.by = 'brainfine')
ggsave("14SingleRvsTabulaM.png", plot = UMAPcelltype1 + UMAPcelltype2  , width = 13, height = 8, dpi = 700)



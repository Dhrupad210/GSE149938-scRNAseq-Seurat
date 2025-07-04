 
# Load Required Libraries
 
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

 
# Load and Subset UMI Matrix
 
umi_matrix <- read.csv("GSE149938_umi_matrix.csv.gz", row.names = 1)

# Subset 3000 random cells to save memory
set.seed(42)
subset_cells <- sample(rownames(umi_matrix), 3000)
subset_matrix <- t(umi_matrix[subset_cells, ])  # Transpose to genes x cells

# Create initial Seurat object
seurat_obj <- CreateSeuratObject(
  counts = subset_matrix,
  project = "GSE149938_sub",
  min.cells = 3,
  min.features = 200
)
saveRDS(seurat_obj, "seurat_step1_raw.rds")

 
# Quality Control
 
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Violin plot for QC metrics
vln <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("vlnplot_qc_metrics.png", plot = vln, width = 10, height = 5)

# Scatter plot: nCount vs nFeature
qc_scatter <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("qc_scatter.png", plot = qc_scatter, width = 6, height = 5)

# Filter cells with too few or too many features
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 250 & nFeature_RNA < 4000)
saveRDS(seurat_obj, "seurat_step2_filtered.rds")

 
#  Normalization and HVG Selection
 
seurat_obj <- NormalizeData(seurat_obj)
saveRDS(seurat_obj, "seurat_step3_normalized.rds")

# Find 2000 highly variable genes (HVGs)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("variable_features_top10.png", plot = plot2, width = 8, height = 5)
saveRDS(seurat_obj, "seurat_step4_variable_features.rds")

 
# Scaling and PCA
 
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
saveRDS(seurat_obj, "seurat_step5_scaled.rds")

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
saveRDS(seurat_obj, "seurat_step6_pca.rds")

# Save PCA visuals
pca_heatmap <- DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
ggsave("pca_heatmap_dim1.png", plot = pca_heatmap, width = 6, height = 5)

elbow <- ElbowPlot(seurat_obj)
ggsave("elbow_plot.png", plot = elbow, width = 6, height = 5)

 
# Clustering and UMAP
 
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
ggsave("step8_umap_clusters.png", plot = umap_plot, width = 8, height = 6)
saveRDS(seurat_obj, "seurat_step8_umap.rds")


# Marker Gene Identification
 
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(all_markers, "step9_all_markers.csv", row.names = FALSE)
saveRDS(all_markers, "step9_all_markers.rds")

 
#  Dot Plot of Top Markers
 
top5 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

dotplot_top5 <- DotPlot(seurat_obj, features = unique(top5$gene)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("dotplot_top5_markers.png", plot = dotplot_top5, width = 12, height = 6)

write.csv(top5, "top5_markers_per_cluster.csv", row.names = FALSE)

 
# Cell Type Annotation
 
celltype_annotations <- c(
  "0" = "Monocytes",      "1" = "B cells",         "2" = "Erythrocytes",
  "3" = "T cells",        "4" = "NK cells",        "5" = "Dendritic",
  "6" = "Progenitors",    "7" = "Plasmablasts",    "8" = "Granulocytes",
  "9" = "Erythroid",      "10" = "pDC",            "11" = "Megakaryocytes",
  "12" = "Naive B",       "13" = "Activated T",    "14" = "Memory B"
)

# Annotate clusters
seurat_obj$celltype <- celltype_annotations[as.character(Idents(seurat_obj))]

# UMAP with cell type labels
celltype_plot <- DimPlot(seurat_obj, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Annotation")
ggsave("umap_celltype_annotation.png", plot = celltype_plot, width = 8, height = 6)

# Save final annotated object
saveRDS(seurat_obj, "seurat_obj_annotated.rds")

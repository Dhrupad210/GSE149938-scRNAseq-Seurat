# 🔬 GSE149938 scRNA-seq Analysis with Seurat

This repository provides a **beginner-friendly pipeline** for analyzing single-cell RNA sequencing (scRNA-seq) data from **GSE149938** using the **Seurat** package in R. The pipeline covers:

- ✅ Raw data loading  
- ✅ Quality control  
- ✅ Normalization  
- ✅ Clustering  
- ✅ Marker gene identification  
- ✅ Cell type annotation  

---

## 📦 Dataset: GSE149938 (Human Blood Cells)

- **Source**: [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149938)  
- **Samples**: 7,643 single-cell transcriptomes from human bone marrow and peripheral blood  
- **Donors**: 21 healthy individuals  
- **Cell Types**: 32 immunophenotypic categories  

---

## 📁 Files Required

| File | Description |
|------|-------------|
| `GSE149938_umi_matrix.csv.gz` | UMI count matrix |
| `GSE149938_barcodes_for_B_NK_Mo_Neu.txt.gz` | Barcodes for B, NK, monocytes, neutrophils |
| `GSE149938_barcodes_for_ery_T.txt.gz` | Barcodes for erythrocytes and T cells |
| `GSM4793029_B.txt.gz`, ..., `GSM4793034_Ery.txt.gz` | Barcodes for individual cell populations |

---

## 🎯 Goal

Identify transcriptionally distinct clusters and **annotate them with known immune cell types** based on marker gene expression.

---

## 🔁 Analysis Workflow

### 🧾 1. Load Raw Data
- Read in the UMI matrix and barcode files
- Handle duplicated gene names
- Check matrix orientation (transpose if cell-by-gene)

### 🔧 2. Create Seurat Object
- Create initial Seurat object
- Calculate mitochondrial % (`percent.mt`)
- **Save**: `seurat_step1_raw.rds`

### 🔍 3. Quality Control & Filtering
- Filter cells with:
  - `nFeature_RNA > 250` and `< 4000`
  - `percent.mt < 5` (optional)
- Visualize QC metrics (violin, scatter)
- **Save**: `seurat_step2_filtered.rds`

### ⚙️ 4. Normalization & Feature Selection
- Normalize using `LogNormalize`
- Identify 2,000 highly variable genes
- **Save**:
  - `seurat_step3_normalized.rds`
  - `seurat_step4_variable_features.rds`

### 📊 5. Scaling & PCA
- Scale entire dataset
- Run PCA on HVGs
- Visualize elbow plot & heatmaps
- **Save**:
  - `seurat_step5_scaled.rds`
  - `seurat_step6_pca.rds`

### 📈 6. Clustering & UMAP
- Find neighbors and clusters (`resolution = 0.5`)
- Run UMAP for 2D visualization
- **Save**:
  - `seurat_step8_umap.rds`
  - `step8_umap_clusters.png`

### 🧬 7. Marker Gene Detection
- Use `FindAllMarkers`:
  - `min.pct = 0.25`, `logfc.threshold = 0.25`
- Identify and export top markers
- **Save**:
  - `step9_all_markers.csv`
  - `top5_markers_per_cluster.csv`
  - `step9_all_markers.rds`

### 📌 8. Cell Type Annotation
- Annotate clusters using canonical marker genes
- Visualize annotated UMAP
- **Save**: `seurat_obj_annotated.rds`

---

## 🎨 Visualization Outputs

| Plot | Description |
|------|-------------|
| `vlnplot_qc_metrics.png` | Violin plot of RNA, gene, and mito content |
| `qc_scatter.png` | Scatter of nCount_RNA vs. nFeature_RNA |
| `elbow_plot.png` | PCA elbow plot |
| `pca_heatmap_dim1.png` | PCA heatmap for first PC |
| `variable_features_top10.png` | Top 10 HVGs |
| `step8_umap_clusters.png` | UMAP with cluster labels |
| `umap_celltype_annotation.png` | UMAP with annotated cell types |
| `dotplot_top5_markers.png` | Dot plot for top markers per cluster |

---

## 📁 Data Outputs

- Intermediate Seurat objects:  
  `seurat_step1_raw.rds`, `seurat_step2_filtered.rds`, ..., `seurat_obj_annotated.rds`

- Marker data:
  - `step9_all_markers.csv`
  - `step9_all_markers.rds`
  - `top5_markers_per_cluster.csv`

---

## 🛠️ How to Run

```
# 1. Clone the Repository
git clone https://github.com/Dhrupad210/GSE149938-scRNAseq-Seurat.git
cd GSE149938-scRNAseq-Seurat

# 2. Download Required Files from GEO (GSE149938)
# Place the following in the working directory:
# - GSE149938_umi_matrix.csv.gz
# - GSE149938_barcodes_for_B_NK_Mo_Neu.txt.gz
# - GSE149938_barcodes_for_ery_T.txt.gz
# - GSM4793029_B.txt.gz to GSM4793034_Ery.txt.gz

# 3. Install Required R Packages
R
install.packages(c("Seurat", "dplyr", "ggplot2", "pheatmap", "Matrix"))
q()

# 4. Run the Analysis Script
Rscript gse149938.R

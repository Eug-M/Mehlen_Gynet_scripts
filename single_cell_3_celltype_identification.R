# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

### from: https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-4-cell-type-identification/

# SingleR + celldex: Reference-based annotation with curated cell type atlases
# scCATCH: Tissue-specific marker databases for automated annotation
# scType dependencies: Tools for gene symbol validation and Excel file reading
# Visualization packages: For creating comparison plots and Sankey diagrams

#-----------------------------------------------
# STEP 1: Load required libraries
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)
library(SeuratObject)

# Automated annotation methods
library(SingleR)
library(celldex)
library(scCATCH)

# Data structures
library(SingleCellExperiment)

# Marker processing
library(HGNChelper)

# Visualization and data manipulation
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(patchwork)
library(scales)
library(viridis)

# Set working directory
setwd("/home/eugenie-modolo/Documents/Gynet/Samples/phase1_singlecell/cell_type_annotation")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/manual_annotation", showWarnings = FALSE)
dir.create("plots/automated_annotation", showWarnings = FALSE)
dir.create("plots/method_comparison", showWarnings = FALSE)
dir.create("annotations", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))



#-----------------------------------------------
# STEP 2: Load integrated data from previous step
#-----------------------------------------------

# Path to integrated Seurat object from previous step
integrated_path <- "/home/eugenie-modolo/Documents/Gynet/Samples/phase1_singlecell/single_cell_integration_analysis/integrated_data/integrated_clustered_seurat.rds"

# Load the object
integrated_seurat <- readRDS(integrated_path)

# Determine which UMAP to use (from best integration method in previous step)
available_reductions <- names(integrated_seurat@reductions)
if ("umap.harmony" %in% available_reductions) {
  umap_reduction <- "umap.harmony"
} else if ("umap.cca" %in% available_reductions) {
  umap_reduction <- "umap.cca"
} else if ("umap.rpca" %in% available_reductions) {
  umap_reduction <- "umap.rpca"
} else if ("umap.mnn" %in% available_reductions) {
  umap_reduction <- "umap.mnn"
} else {
  umap_reduction <- "umap"
}



#-----------------------------------------------
# STEP 3: Visualize clusters before annotation
#-----------------------------------------------

# Create multi-panel overview
p1 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "seurat_clusters", label = TRUE, label.size = 5,
              pt.size = 0.1) +
  ggtitle("Clusters (Pre-Annotation)") +
  theme(plot.title = element_text(face = "bold", size = 14))

p2 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "sample_id", pt.size = 0.1) +
  ggtitle("Samples") +
  theme(legend.text = element_text(size = 7))

p3 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "condition", pt.size = 0.1) +
  ggtitle("Condition") +
  scale_color_manual(values = c("Healthy" = "#2E86AB", 
                                "Post_Treatment" = "#F18F01"))

# Combine
p_overview <- (p1 | p2) / (p3 | plot_spacer())

ggsave("plots/00_starting_clusters.png", p_overview, 
       width = 14, height = 10, dpi = 300)



#-----------------------------------------------
# STEP 4: Normalize and scale data (Critical: FindAllMarkers, FeaturePlot, DotPlot, and DoHeatmap all require normalized data)
#-----------------------------------------------

# Normalize data (log-normalization)
integrated_seurat <- NormalizeData(integrated_seurat, 
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000, 
                                   verbose = FALSE)

# Scale all genes (required for DoHeatmap and marker analysis)
integrated_seurat <- ScaleData(integrated_seurat, 
                               features = rownames(integrated_seurat), 
                               verbose = FALSE)



#-----------------------------------------------
# STEP 5: Define and check canonical markers
#-----------------------------------------------

# Define comprehensive marker panel for PBMC cell types
canonical_markers <- list(
  "immune_cells" = c("PTPRC", "CD86", "CD4", "CD8A", "CD28", "CD68"),
  "lymphocytes" = c("CD247", "CD3G", "CD3D", "CD3E", "CD2"),
  "dendrits_and_monocytes" = c("ITGAX", "CD14", "MAFB"),
  "Treg_lymphocytes" = c("CD2", "FOXP3", "CDH5", "CD8A"),
  "activated_lymphocytes" = c("GZMB", "GZMA"),
  "granzymes" = c("FASLG", "TGFBR1", "ITGB8", "CD4"),
  "Epithelial" = c("EPCAM"),
  "Pericyte" = c("ACTA2", "CSPG4"),
  "Endothelium" = c("CDH5", "PECAM1"),
  "Endothelium_lymphatic" = c("PROX1", "LYVE1"),
  # "" = c(),
  "T_cells" = c("CD3D", "CD3E", "CD3G"),
  "CD4_T" = c("CD3D", "CD4", "IL7R"),
  "CD8_T" = c("CD3D", "CD8A", "CD8B"),
  "B_cells" = c("CD79A", "MS4A1", "CD19"),
  "NK_cells" = c("NKG7", "GNLY", "NCAM1"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Classical_Mono" = c("CD14", "S100A8", "FCGR3A"),
  "Non_Classical_Mono" = c("FCGR3A", "MS4A7"),
  # "Dendritic_cells" = c("FCER1A", "CST3"),
  # "Plasmacytoid_DC" = c("IL3RA", "GZMB", "SERPINF1"),
  "Platelets" = c("PPBP", "PF4", "GP9")
)

# Check which markers are present in dataset
all_genes <- rownames(integrated_seurat)

for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  present <- markers %in% all_genes
  
  cat(sprintf("%-20s: %d/%d markers present\n", 
              cell_type, sum(present), length(markers)))
  
  if (!all(present)) {
    missing <- markers[!present]
    cat(sprintf("  Missing: %s\n", paste(missing, collapse = ", ")))
  }
}



#-----------------------------------------------
# STEP 6: Visualize canonical markers
#-----------------------------------------------

# Function to create violin plots for a marker set
plot_violin_markers <- function(markers, cell_type_name, filename_suffix) {
  present_markers <- markers[markers %in% rownames(integrated_seurat)]
  
  if (length(present_markers) == 0) {
    cat("  ⚠ No markers present for", cell_type_name, "\n")
    return(NULL)
  }
  
  cat("  •", cell_type_name, ":", length(present_markers), "markers\n")
  
  p <- VlnPlot(integrated_seurat,
               features = present_markers,
               group.by = "seurat_clusters",
               pt.size = 0,
               ncol = 3) &
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold"))
  
  n_rows <- ceiling(length(present_markers) / 3)
  height <- max(4, n_rows * 2.5)
  
  ggsave(paste0("plots/manual_annotation/violin_", filename_suffix, ".png"),
         p, width = 14, height = height, dpi = 300)
  
  return(p)
}

# Function to create UMAP feature plots for a marker set
plot_umap_markers <- function(markers, cell_type_name, filename_suffix) {
  present_markers <- markers[markers %in% rownames(integrated_seurat)]
  
  if (length(present_markers) == 0) {
    return(NULL)
  }
  
  p <- FeaturePlot(integrated_seurat,
                   features = present_markers,
                   reduction = umap_reduction,
                   ncol = 4,
                   pt.size = 0.1) &
    theme(plot.title = element_text(size = 10),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right")
  
  n_rows <- ceiling(length(present_markers) / 4)
  height <- max(4, n_rows * 3)
  
  ggsave(paste0("plots/manual_annotation/umap_", filename_suffix, ".png"),
         p, width = 16, height = height, dpi = 300)
  
  return(p)
}

# Generate plots for ALL cell types
for (cell_type_key in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type_key]]
  cell_type_name <- gsub("_", " ", cell_type_key)
  
  plot_violin_markers(markers, cell_type_name, cell_type_key)
  plot_umap_markers(markers, cell_type_name, cell_type_key)
}

# Create comprehensive DotPlot showing all key markers
all_markers <- unique(unlist(canonical_markers))
all_markers_present <- all_markers[all_markers %in% rownames(integrated_seurat)]

p_dotplot <- DotPlot(integrated_seurat,
                     features = all_markers_present,
                     group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8)) +
  labs(title = "All Canonical Markers by Cluster",
       subtitle = paste(length(all_markers_present), "markers across all cell types"))

ggsave("plots/manual_annotation/dotplot_all_canonical_markers.png", p_dotplot,
       width = 12, height = max(8, length(all_markers_present) * 0.15), dpi = 300)



#-----------------------------------------------
# STEP 7: Find cluster-specific marker genes
#-----------------------------------------------

# Find markers for each cluster vs all other cells
cluster_markers <- FindAllMarkers(
  integrated_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

# Remove ribosomal and mitochondrial genes
cluster_markers_filtered <- cluster_markers %>%
  filter(!grepl("^RP[SL]", gene)) %>%
  filter(!grepl("^MT-", gene))

# Get top 5 markers per cluster
top_markers <- cluster_markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

cat("\nTop 5 markers per cluster:\n")
print(top_markers %>% select(cluster, gene, avg_log2FC, pct.1, pct.2))

# Save all markers
write.csv(cluster_markers_filtered, 
          "annotations/cluster_markers_all.csv",
          row.names = FALSE)

write.csv(top_markers,
          "annotations/top5_markers_per_cluster.csv",
          row.names = FALSE)

# Create heatmap of top markers
top_genes <- top_markers$gene

p_heatmap <- DoHeatmap(
  integrated_seurat,
  features = top_genes,
  group.by = "seurat_clusters",
  size = 3
) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "Top 5 Markers per Cluster")

ggsave("plots/manual_annotation/06_heatmap_top_markers.png", p_heatmap,
       width = 12, height = 14, dpi = 300)



#-----------------------------------------------
# STEP 8: Manual cell type assignment
#-----------------------------------------------

# /!\ à ajouter : étape inferCNV pour les clusters tumoraux

# MAPPING - CUSTOMIZED using the figures from Step 7
cluster_to_celltype <- c(
  "0" = "lymphocytes",
  "1" = "monocytes",
  "2" = "epithelial cells",
  "3" = "epithelial cells",
  "4" = "blood cells",
  "5" = "CD8+ T cells",
  "6" = "monocytes",
  "7" = "CD4+ T cells",
  "8" = "NK cells",
  "9" = "endothelial cells",
  "10" = "pericyte cells",
  "11" = "CD4+ T cells",
  "12" = "epithelial_tumor",
  "13" = "CD4+ T cells",
  "14" = "epithelial_tumor",
  "15" = "B Cells",
  "16" = "epithelial_tumor",
  "17" = "epithelial cells"
)

# Apply to all cells
cluster_ids <- as.character(integrated_seurat$seurat_clusters)
integrated_seurat$manual_annotation <- unname(cluster_to_celltype[cluster_ids])

# Visualize manual annotation
p_manual <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("Manual Cell Type Annotation") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/manual_annotation/07_manual_annotation_umap.png", p_manual,
       width = 10, height = 8, dpi = 300)

# Split by condition
p_manual_split <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  split.by = "condition",
  pt.size = 0.1,
  ncol = 2
) +
  ggtitle("Manual Annotation by Condition")

ggsave("plots/manual_annotation/08_manual_annotation_by_condition.png",
       p_manual_split, width = 14, height = 6, dpi = 300)



#-----------------------------------------------
# STEP 9: SingleR with multiple reference datasets
#-----------------------------------------------

# Convert Seurat to SingleCellExperiment for SingleR
sce <- as.SingleCellExperiment(integrated_seurat)

# Reference 1: HumanPrimaryCellAtlas (Broad)
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

singler_hpca <- SingleR(
  test = sce,
  ref = hpca_ref,
  labels = hpca_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_hpca <- singler_hpca$labels

# Reference 2: MonacoImmuneData (Immune-Specific)
monaco_ref <- celldex::MonacoImmuneData()

singler_monaco <- SingleR(
  test = sce,
  ref = monaco_ref,
  labels = monaco_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_monaco <- singler_monaco$labels

# Reference 3: DatabaseImmuneCellExpression
dice_ref <- celldex::DatabaseImmuneCellExpressionData()

singler_dice <- SingleR(
  test = sce,
  ref = dice_ref,
  labels = dice_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_dice <- singler_dice$labels

# Reference 4: BlueprintEncodeData
dice_ref <- celldex::BlueprintEncodeData()

singler_bed <- SingleR(
  test = sce,
  ref = dice_ref,
  labels = dice_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_bed <- singler_bed$labels



#-----------------------------------------------
# STEP 10: Compare SingleR references
#-----------------------------------------------

# Visualize all three references
p_hpca <- DimPlot(integrated_seurat, reduction = umap_reduction,
                  group.by = "singler_hpca", pt.size = 0.1) +
  labs(title = "HPCA Reference") +
  theme(legend.text = element_text(size = 7))

p_monaco <- DimPlot(integrated_seurat, reduction = umap_reduction,
                    group.by = "singler_monaco", pt.size = 0.1) +
  labs(title = "Monaco Reference") +
  theme(legend.text = element_text(size = 7))
# 
# p_dice <- DimPlot(integrated_seurat, reduction = umap_reduction,
#                   group.by = "singler_dice", pt.size = 0.1) +
#   labs(title = "DICE Reference") +
#   theme(legend.text = element_text(size = 7))

p_bed <- DimPlot(integrated_seurat, reduction = umap_reduction,
                  group.by = "singler_bed", pt.size = 0.1) +
  labs(title = "BED Reference") +
  theme(legend.text = element_text(size = 7))

p_manual_ref <- DimPlot(integrated_seurat, reduction = umap_reduction,
                        group.by = "manual_annotation", pt.size = 0.1) +
  labs(title = "Manual") +
  theme(legend.text = element_text(size = 7))

p_ref_comparison <- (p_manual_ref | p_hpca) / (p_monaco | p_bed)

ggsave("plots/automated_annotation/10_singler_reference_comparison.png",
       p_ref_comparison, width = 16, height = 12, dpi = 300)

# Choosing the annotation
integrated_seurat$singler_annotation <- integrated_seurat$singler_bed

# Visualize chosen reference
p_singler_final <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "singler_annotation",
  label = TRUE,
  label.size = 3,
  pt.size = 0.1,
  repel = TRUE
) +
  ggtitle("SingleR Annotation (Monaco Reference)") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/automated_annotation/11_singler_final_annotation.png",
       p_singler_final, width = 11, height = 8, dpi = 300)



# /!\ continue: from step 11




### Code Nico pour trouver les cellules tumorales, une fois les cellules immunitaires trouvées 
#https://github.com/broadinstitute/infercnv/wiki/infercnv-10x
library(infercnv)
rawcount <- as.matrix(LayerData(obj, assay = "RNA", layer = "counts"))

#sampleannotations
sample.annotation <- obj@meta.data[, c("celltype", "orig.ident")]
sample.annotation$orig.ident <- rownames(sample.annotation)
head(sample.annotation)
colnames(sample.annotation) <- NULL


geneorder <- read.table("/Users/nicolasrama/Desktop/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions.txt", header = T)
head(geneorder)
rownames(geneorder) <- geneorder$GeneId
geneorder <- geneorder[, -1]
head(geneorder)

table(colnames(rawcount) == rownames(sample.annotation))

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = rawcount, annotations_file = sample.annotation, gene_order_file = geneorder, ref_group_names = c("Neutrophils", "Mono/macrophage", "lymphocyteT", "LymphocyteB"))



infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(),
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE, num_threads = 1)


seurat_obj = infercnv::add_to_seurat(infercnv_output_path='/var/folders/j2/0rkf_l2x10z866tvy9bss1_h0000gn/T/RtmpB4z3ET/file58e1e24729f/',
                                     seurat_obj=obj, # optional
                                     top_n=10
)


FeaturePlot(seurat_obj, features = "proportion_scaled_cnv_chr19", order = T, label=F, split.by = "orig.ident") #
FeaturePlot(seurat_obj, features = "mitoRatio", order = T, label=T, min.cutoff = "q5") #

Idents(obj) <- obj@meta.data$SCT_snn_res.0.2
DimPlot(obj, split.by = "orig.ident")
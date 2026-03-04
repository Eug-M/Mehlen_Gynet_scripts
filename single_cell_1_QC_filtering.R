# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

### from: https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-2-quality-control-and-cell-filtering/

#-----------------------------------------------
# STEP 1: Load libraries and configure environment
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)              # Main scRNA-seq toolkit
library(SeuratObject)        # Seurat object infrastructure

# QC-specific packages
library(DropletUtils)        # Empty droplet detection
library(scDblFinder)         # Doublet detection
library(SoupX)               # Ambient RNA correction
library(SingleCellExperiment) # SCE objects

# Visualization
library(ggplot2)             # Custom plots when needed
library(patchwork)           # Combine plots
library(dplyr)               # Data manipulation

# Set working directory
path <- "/home/eugenie-modolo/Documents/Gynet/Samples/phase1_singlecell/single_cell_QC_analysis"
setwd(path)

# # Create output directories
# dir.create("plots", showWarnings = FALSE)
# dir.create("qc_metrics", showWarnings = FALSE)
# dir.create("filtered_data", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Verify Seurat 5 installation
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

# Check for Seurat 5 features
if (packageVersion("Seurat") >= "5.0.0") {
  cat("✓ Seurat 5 detected - layer-based architecture available\n")
} else {
  warning("Please upgrade to Seurat 5 for this tutorial")
}


#-----------------------------------------------
# STEP 2: Load 10x Genomics data
#-----------------------------------------------

# Define paths
# sample_name <- "C1D1"
sample_name <- "C3D1"
cellranger_output <- "../"

# OPTION A: Load RAW matrix (recommanded)
# counts <- Read10X(data.dir = file.path(cellranger_output, "raw_feature_bc_matrix"))

# OPTION B: Load FILTERED matrix (alternative - skip Steps 4-5 if using this)
counts <- Read10X(data.dir = file.path(cellranger_output, paste0(sample_name, "_cellranger")))

# Create Seurat object
# If raw: contains ALL droplets (~20,000-100,000)
# If filtered: contains only cells Cell Ranger called (~3,000-5,000)
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = sample_name,
  min.cells = 0,      # Don't filter genes yet
  min.features = 0    # Don't filter cells yet - we'll do this properly
)

# Add sample metadata (consistent with Part 1 tutorial)
seurat_obj$sample_id <- sample_name
seurat_obj$sra_id <- "SRRxxx"
seurat_obj$condition <- sample_name
seurat_obj$patient_id <- paste0("Donor_", sample_name)
seurat_obj$time_point <- NA  # Not applicable for healthy donors

# Quick check
cat("Loaded:", ncol(seurat_obj), "droplets ×", nrow(seurat_obj), "genes\n")
cat("(Most are empty - EmptyDrops will filter in Step 4)\n")


#-----------------------------------------------
# STEP 3: Initial data exploration
#-----------------------------------------------

# Basic statistics
cat("Total droplets:", ncol(seurat_obj), "\n")
cat("Total genes:", nrow(seurat_obj), "\n")

# Quick statistics
sparsity <- 1 - (sum(LayerData(seurat_obj, layer = "counts") > 0) / 
                   (nrow(seurat_obj) * ncol(seurat_obj)))
cat("Matrix sparsity:", round(sparsity * 100, 1), "%\n")

# UMI distribution - will be bimodal (empty droplets + real cells)
umi_counts <- colSums(LayerData(seurat_obj, layer = "counts"))
cat("Median UMI per droplet:", median(umi_counts), "\n")
cat("Droplets with >500 UMI:", sum(umi_counts > 500), 
    "(likely cells, will be validated by EmptyDrops)\n")


# #-----------------------------------------------
# # STEP 4: Empty droplet detection - skipped for GyNet (we don't have the raw data before Cell Ranger)
# #-----------------------------------------------
# 
# Convert Seurat object to SingleCellExperiment for EmptyDrops
# We already loaded the raw matrix in Step 2
sce <- as.SingleCellExperiment(seurat_obj)

# Track initial droplet count for QC summary
initial_droplet_count <- ncol(sce)
initial_gene_count <- nrow(sce)
# 
# # Run EmptyDrops
# set.seed(100)
# empty_results <- emptyDrops(
#   m = counts(sce),
#   lower = 100,        # Droplets with <100 UMI assumed empty (for estimating ambient RNA)
#   niters = 10000,     # Iterations for p-value calculation
#   test.ambient = TRUE # Test if droplets differ from ambient RNA profile
# )
# 
# # Identify cells (FDR < 0.01)
# is_cell <- empty_results$FDR < 0.01
# is_cell[is.na(is_cell)] <- FALSE  # Treat NA as empty (very low UMI droplets)
# validated_barcodes <- colnames(sce)[is_cell]
# 
# # Brief feedback
# cat("Droplets tested:", ncol(sce), "\n")
# cat("Cells called:", sum(is_cell), 
#     "(", round(sum(is_cell)/ncol(sce)*100, 1), "%)\n")
# cat("Empty droplets removed:", sum(!is_cell), 
#     "(", round(sum(!is_cell)/ncol(sce)*100, 1), "%)\n")
# 
# # Visual QC - the critical check!
# empty_df <- data.frame(
#   total_umi = colSums(counts(sce)),
#   is_cell = is_cell
# )
# 
# p1 <- ggplot(empty_df, aes(x = log10(total_umi + 1), fill = is_cell)) +
#   geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
#   scale_fill_manual(
#     values = c("TRUE" = "#2E86AB", "FALSE" = "#A23B72"),
#     labels = c("Empty Droplet", "Cell"),
#     name = "Classification"
#   ) +
#   labs(
#     title = "EmptyDrops: Cell vs Empty Droplet Detection",
#     subtitle = paste(sum(is_cell), "cells called from", ncol(sce), "total droplets"),
#     x = "log10(UMI + 1)",
#     y = "Number of Droplets"
#   ) +
#   theme_classic() +
#   theme(plot.title = element_text(face = "bold", size = 14))
# 
# ggsave("plots/01_empty_droplets.png", p1, width = 10, height = 6, dpi = 300)
# 
# # NOW filter the Seurat object to keep only validated cells
# # IMPORTANT: Save raw counts BEFORE filtering (needed for SoupX in Step 5)
# raw_counts_all_droplets <- LayerData(seurat_obj, layer = "counts")
# 
# seurat_obj <- subset(seurat_obj, cells = validated_barcodes)
# 
# cat("After EmptyDrops:", ncol(seurat_obj), "cells retained\n")
# 
# 
# #-----------------------------------------------
# # STEP 5: Ambient RNA correction with SoupX - skipped for GyNet (we don't have the raw data before Cell Ranger)
# #-----------------------------------------------
# 
# # Prepare SoupX channel
# # tod (table of droplets) = raw counts from ALL droplets (saved in Step 4)
# # toc (table of cells) = filtered counts from validated cells only
# tod <- raw_counts_all_droplets
# toc <- LayerData(seurat_obj, layer = "counts")
# sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)
# 
# # Quick clustering for contamination estimation
# temp_obj <- seurat_obj
# temp_obj <- NormalizeData(temp_obj, verbose = FALSE)
# temp_obj <- FindVariableFeatures(temp_obj, nfeatures = 2000, verbose = FALSE)
# temp_obj <- ScaleData(temp_obj, verbose = FALSE)
# temp_obj <- RunPCA(temp_obj, npcs = 30, verbose = FALSE)
# temp_obj <- FindNeighbors(temp_obj, dims = 1:30, verbose = FALSE)
# temp_obj <- FindClusters(temp_obj, resolution = 0.8, verbose = FALSE)
# sc <- setClusters(sc, setNames(as.character(temp_obj$seurat_clusters), colnames(temp_obj)))
# 
# # Estimate contamination
# sc <- tryCatch({
#   autoEstCont(sc, verbose = FALSE)
# }, error = function(e) {
#   cat("Note: autoEstCont failed\n")
#   sc$fit$rho <- NULL
#   return(sc)
# })
# 
# contamination_fraction <- sc$fit$rho
# 
# # Check contamination level
# cat("Contamination fraction:", 
#     ifelse(is.null(contamination_fraction), "NULL (very clean sample)", 
#            paste0(round(contamination_fraction * 100, 2), "%")), "\n")
# 
# # Decision: Skip or correct based on contamination level
# if (is.null(contamination_fraction) || contamination_fraction < 0.05) {
#   cat("→ Contamination is NULL or <5% - skipping SoupX correction\n")
#   cat("   Your sample is clean! Proceeding with original counts.\n")
# } else {
#   cat("→ Contamination is", round(contamination_fraction * 100, 2), "% (≥5%)\n")
#   cat("   Applying SoupX correction...\n")
#   
#   # Apply correction
#   suppressWarnings({
#     corrected_counts <- adjustCounts(sc)
#   })
#   seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = corrected_counts)
#   cat("   ✓ SoupX correction applied\n")
# }



#-----------------------------------------------
# STEP 6: Doublet detection with scDblFinder
#-----------------------------------------------

# Join layers and convert to SCE
seurat_obj <- JoinLayers(seurat_obj)
sce <- as.SingleCellExperiment(seurat_obj)

# Run scDblFinder (suppress expected warnings about empty layers and xgboost)
set.seed(42)
suppressWarnings({
  sce <- scDblFinder(sce, dbr = NULL, verbose = FALSE)
})

# Extract results
seurat_obj$doublet_class <- ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet")
seurat_obj$doublet_score <- sce$scDblFinder.score

# Brief feedback
n_doublets <- sum(seurat_obj$doublet_class == "Doublet")
doublet_rate <- n_doublets / ncol(seurat_obj) * 100
cat("Doublets:", n_doublets, "(", round(doublet_rate, 1), "%)\n")

# Track cells before doublet removal for QC summary
cells_before_doublet_removal <- ncol(seurat_obj)

# Quick UMAP for visualization
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

# Visualize doublets on UMAP
p3 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "doublet_class",
  cols = c("Singlet" = "#90E0EF", "Doublet" = "#EF233C"),
  pt.size = 1
) +
  labs(title = paste0("Doublet Detection (", round(doublet_rate, 1), "%)")) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(paste0("plots/", sample_name, "_02_doublets_umap.png"), p3, width = 8, height = 7, dpi = 300)

# Remove doublets
seurat_obj <- subset(seurat_obj, subset = doublet_class == "Singlet")
cat("After doublet removal:", ncol(seurat_obj), "cells\n")



#-----------------------------------------------
# STEP 7: Calculate QC metrics, visualize, and filter cells
#-----------------------------------------------

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# Summary statistics to guide threshold setting
cat("\nQC Metric Distributions:\n")
cat("nCount_RNA (UMI):\n")
cat("  Median:", median(seurat_obj$nCount_RNA), 
    "| Q1-Q3:", quantile(seurat_obj$nCount_RNA, 0.25), "-", 
    quantile(seurat_obj$nCount_RNA, 0.75), "\n")

cat("nFeature_RNA (genes):\n")
cat("  Median:", median(seurat_obj$nFeature_RNA),
    "| Q1-Q3:", quantile(seurat_obj$nFeature_RNA, 0.25), "-",
    quantile(seurat_obj$nFeature_RNA, 0.75), "\n")

cat("percent.mt:\n")
cat("  Median:", round(median(seurat_obj$percent.mt), 2), "%",
    "| 95th percentile:", round(quantile(seurat_obj$percent.mt, 0.95), 2), "%\n")

cat("percent.ribo:\n")
cat("  Median:", round(median(seurat_obj$percent.ribo), 2), "%",
    "| 95th percentile:", round(quantile(seurat_obj$percent.ribo, 0.95), 2), "%\n")

# Visualize distributions - THE CRITICAL STEP
p4 <- VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  pt.size = 0.1 # 0 to remove the dots, 0.1 to see them
) &
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0("plots/", sample_name, "_03_qc_violins.png"), p4, width = 16, height = 4, dpi = 300)

# Scatter plots showing metric relationships
p5 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  labs(title = "UMI vs Genes Detected")

p6 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  labs(title = "UMI vs Mitochondrial %")

p7 <- FeatureScatter(seurat_obj, feature1 = "percent.mt", feature2 = "percent.ribo") +
  labs(title = "Mitochondrial % vs Ribosomal %")

p_scatter <- p5 + p6 + p7
ggsave(paste0("plots/", sample_name, "_04_qc_scatter.png"), p_scatter, width = 15, height = 5, dpi = 300)

cat("\n→ EXAMINE plots/03_qc_violins.png and plots/04_qc_scatter.png\n")
cat("→ Look for:\n")
cat("   - Bimodal distributions (good vs bad cells)\n")
cat("   - Outliers in scatter plots\n")
cat("   - Relationship between metrics\n\n")

# Set thresholds based on visual inspection
# MANUALLY ADJUST THESE based on what you see in the plots!

cat("=== Setting Filtering Thresholds ===\n")
cat("Based on visual inspection of the plots above:\n\n")

# For our case - adjust these for your specific sample!
# C1D1
# nfeature_min <- 300    # Minimum genes per cell (aussi testé 400)
# nfeature_max <- 7000   # Maximum genes (above this = likely doublets) 
# ncount_min <- 600      # Minimum UMI per cell (aussi testé 700)
# ncount_max <- 30000    # Maximum UMI (above this = likely doublets) (aussi testé 25000)
# mt_thresh <- 25        # Maximum mitochondrial % (adjust based on tissue) (aussi testé 30)
# C3D1
nfeature_min <- 300    # Minimum genes per cell
nfeature_max <- 5000   # Maximum genes (above this = likely doublets)
ncount_min <- 500      # Minimum UMI per cell
ncount_max <- 30000    # Maximum UMI (above this = likely doublets)
mt_thresh <- 25        # Maximum mitochondrial % (adjust based on tissue)

cat("Set thresholds (adjust based on YOUR data):\n")
cat("  nFeature_RNA: [", nfeature_min, ",", nfeature_max, "]\n")
cat("  nCount_RNA: [", ncount_min, ",", ncount_max, "]\n")
cat("  percent.mt: <", mt_thresh, "%\n\n")

cat("Guidelines for threshold setting:\n")
cat("  • nFeature_RNA min: Where does the low-quality tail end? (typically 500-1000)\n")
cat("  • nFeature_RNA max: Where do potential doublets start? (typically 4000-7000)\n")
cat("  • percent.mt: PBMCs <10%, Brain <5%, Tumor <20%\n")
cat("  • If distributions are bimodal, set threshold between modes\n")
cat("  • Aim to remove no more than 25% of cells \n\n")

# Visualize thresholds on data
qc_df <- data.frame(
  nCount_RNA = seurat_obj$nCount_RNA,
  nFeature_RNA = seurat_obj$nFeature_RNA,
  percent.mt = seurat_obj$percent.mt
)

qc_df$pass_qc <- (
  qc_df$nCount_RNA >= ncount_min &
    qc_df$nCount_RNA <= ncount_max &
    qc_df$nFeature_RNA >= nfeature_min &
    qc_df$nFeature_RNA <= nfeature_max &
    qc_df$percent.mt < mt_thresh
)

p8 <- ggplot(qc_df, aes(x = log10(nCount_RNA + 1), y = log10(nFeature_RNA + 1), 
                        color = pass_qc)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = log10(c(ncount_min, ncount_max)), 
             linetype = "dashed", color = "red") +
  geom_hline(yintercept = log10(c(nfeature_min, nfeature_max)), 
             linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#06D6A0", "FALSE" = "#EF476F")) +
  labs(
    title = "Cell Filtering Thresholds",
    subtitle = paste0(sum(qc_df$pass_qc), " cells pass QC (", 
                      round(sum(qc_df$pass_qc)/nrow(qc_df)*100, 1), "%)"),
    x = "log10(UMI + 1)",
    y = "log10(Genes + 1)",
    color = "Pass QC"
  ) +
  theme_classic()

ggsave(paste0("plots/", sample_name, "_05_filtering_thresholds.png"), p8, width = 8, height = 7, dpi = 300)

# Report filtering impact
removal_pct <- sum(!qc_df$pass_qc) / nrow(qc_df) * 100
cat("Filtering impact:\n")
cat("  Cells before:", nrow(qc_df), "\n")
cat("  Cells passing QC:", sum(qc_df$pass_qc), "\n")
cat("  Cells removed:", sum(!qc_df$pass_qc), 
    "(", round(removal_pct, 1), "%)\n\n")

# Breakdown by criterion
cat("Removal breakdown:\n")
cat("  Low genes (<", nfeature_min, "):", 
    sum(qc_df$nFeature_RNA < nfeature_min), "\n")
cat("  High genes (>", nfeature_max, "):", 
    sum(qc_df$nFeature_RNA > nfeature_max), "\n")
cat("  High MT% (>", mt_thresh, "%):", 
    sum(qc_df$percent.mt > mt_thresh), "\n\n")
# Check the breakdown – if >90% of removals are from one criterion, that threshold is wrong

# Sanity check on removal rate
if (removal_pct > 30) {
  cat("⚠ WARNING: Removing", round(removal_pct, 1), "% of cells is high!\n")
  cat("   Consider relaxing thresholds (especially MT%)\n\n")
} else if (removal_pct < 5) {
  cat("⚠ WARNING: Only removing", round(removal_pct, 1), "% of cells is low.\n")
  cat("   Check if you're retaining low-quality cells.\n\n")
} else {
  cat("✓ Removal rate is reasonable (typical: 10-25%)\n\n")
}

# Track cells before cell-level QC filtering for summary
cells_before_cell_qc <- nrow(qc_df)

# Apply filters
seurat_obj <- subset(
  seurat_obj,
  subset = nCount_RNA >= ncount_min &
    nCount_RNA <= ncount_max &
    nFeature_RNA >= nfeature_min &
    nFeature_RNA <= nfeature_max &
    percent.mt < mt_thresh
)

cat("After filtering:", ncol(seurat_obj), "cells remaining\n")



#-----------------------------------------------
# STEP 8: Gene-level QC with detection threshold
#-----------------------------------------------

# Calculate gene detection
counts_matrix <- LayerData(seurat_obj, layer = "counts")
gene_detection <- rowSums(counts_matrix > 0)

gene_qc <- data.frame(
  gene = rownames(seurat_obj),
  n_cells_detected = gene_detection,
  pct_cells_detected = (gene_detection / ncol(seurat_obj)) * 100,
  is_mt = grepl("^MT-", rownames(seurat_obj)),
  is_ribo = grepl("^RP[SL]", rownames(seurat_obj)),
  is_hb = grepl("^HB[AB]", rownames(seurat_obj))
)

cat("Total genes:", nrow(gene_qc), "\n")
cat("MT genes:", sum(gene_qc$is_mt), 
    "| Ribo genes:", sum(gene_qc$is_ribo), 
    "| Hb genes:", sum(gene_qc$is_hb), "\n\n")

# Visualize gene detection
p9 <- ggplot(gene_qc, aes(x = pct_cells_detected)) +
  geom_histogram(bins = 50, fill = "#118AB2", alpha = 0.7) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(
    title = "Gene Detection Across Cells",
    subtitle = "Red line: 0.1% of cells threshold",
    x = "% of Cells Expressing Gene (log scale)",
    y = "Number of Genes"
  ) +
  theme_classic()

ggsave(paste0("plots/", sample_name, "_06_gene_detection.png"), p9, width = 8, height = 5, dpi = 300)

cat("→ EXAMINE plots/06_gene_detection.png\n")
cat("   Choose threshold based on where detection drops off\n\n")

# Set filtering threshold - ADJUST BASED ON YOUR DATA
# Common approaches:
#   - 0.1% of cells (very lenient, keeps most genes)
#   - 1% of cells (moderate, standard approach)
#   - 3 cells minimum (absolute count, conservative)

min_pct_cells <- 0.1  # Gene must be in ≥0.1% of cells
# Alternatively: min_cells <- 3  # Gene must be in ≥3 cells

cat("Setting gene filter threshold:\n")
cat("  Minimum detection: ≥", min_pct_cells, "% of cells\n")
cat("  (Equivalent to ≥", ceiling(ncol(seurat_obj) * min_pct_cells / 100), 
    "cells with current dataset)\n\n")

cat("Threshold options:\n")
cat("  • 0.1% of cells: Lenient (keeps rare cell type markers)\n")
cat("  • 1% of cells: Standard (balances detection vs noise)\n")
cat("  • 3-10 cells minimum: Conservative (removes very rare genes)\n\n")

# Filter genes
genes_to_keep <- (gene_qc$pct_cells_detected >= min_pct_cells) & !gene_qc$is_hb

cat("Genes passing filter:", sum(genes_to_keep), "/", nrow(gene_qc), "\n")
cat("Genes removed:\n")
cat("  Low detection:", sum(gene_qc$pct_cells_detected < min_pct_cells), "\n")
cat("  Hemoglobin:", sum(gene_qc$is_hb), "\n\n")

seurat_obj <- seurat_obj[gene_qc$gene[genes_to_keep], ]

cat("After gene filtering:", nrow(seurat_obj), "genes remaining\n")

p10 <- VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  pt.size = 0.1 # 0 to remove the dots
) &
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0("plots/", sample_name, "_03_qc_violins_after.png"), p10, width = 16, height = 4, dpi = 300)

# Scatter plots showing metric relationships
p11 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  labs(title = "UMI vs Genes Detected")

p12 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  labs(title = "UMI vs Mitochondrial %")

p13 <- FeatureScatter(seurat_obj, feature1 = "percent.mt", feature2 = "percent.ribo") +
  labs(title = "Mitochondrial % vs Ribosomal %")

p_scatter <- p11 + p12 + p13
ggsave(paste0("plots/", sample_name, "_04_qc_scatter_after.png"), p_scatter, width = 15, height = 5, dpi = 300)



#-----------------------------------------------
# STEP 9: Normalization and variable features
#-----------------------------------------------

# Log-normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                            scale.factor = 10000, verbose = FALSE)

# Find highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)

top10 <- head(VariableFeatures(seurat_obj), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")
# Top 10 variable genes: 
# C1D1 : 
# C3D1 : TPSB2, IGHG3, IGLV2-11, IGHV1-18, IGLC1, CST1, CCL18, CST4, IGHV3-7, IGLC3

# Visualize variable features
p14 <- VariableFeaturePlot(seurat_obj)
p14 <- LabelPoints(plot = p14, points = top10, repel = TRUE)
ggsave(paste0("plots/", sample_name, "_07_variable_features.png"), p14, width = 10, height = 7, dpi = 300)



#-----------------------------------------------
# STEP 10: Save processed data and create summary
#-----------------------------------------------

# Save the clean Seurat object
saveRDS(seurat_obj, file = paste0("filtered_data/", sample_name, "_qc_filtered.rds"))

# Create ONE comprehensive QC summary (instead of many CSV files)
qc_summary <- data.frame(
  sample = sample_name,
  
  # Cell counts through QC pipeline
  initial_droplets = initial_droplet_count,  # From raw matrix before EmptyDrops
  # after_emptydrops = length(validated_barcodes),  # Cells called by EmptyDrops, Step 4
  # after_soupx = length(validated_barcodes),       # SoupX doesn't remove cells, Step 4
  after_doublets = cells_before_doublet_removal - n_doublets,  # After doublet removal
  after_cell_qc = ncol(seurat_obj),  # After cell-level QC filtering
  final_cells = ncol(seurat_obj),    # Final count (same as after_cell_qc)
  
  # Gene counts
  initial_genes = initial_gene_count,  # From raw matrix
  final_genes = nrow(seurat_obj),      # After gene filtering
  
  # Key metrics (from final filtered data)
  median_umi = median(seurat_obj$nCount_RNA),
  median_genes = median(seurat_obj$nFeature_RNA),
  median_mt_pct = median(seurat_obj$percent.mt),
  
  # QC parameters
  # contamination_pct = ifelse(is.null(contamination_fraction), 0, 
  #                            round(contamination_fraction * 100, 2)),  # Step 5
  doublet_rate_pct = round(doublet_rate, 2),
  
  # Filtering thresholds used
  ncount_min = ncount_min,
  ncount_max = ncount_max,
  nfeature_min = nfeature_min,
  nfeature_max = nfeature_max,
  mt_threshold = mt_thresh
)

# Save comprehensive summary
write.csv(qc_summary, paste0("qc_metrics/", sample_name, "_QC_summary.csv"), row.names = FALSE)
write.csv(seurat_obj@meta.data, paste0("filtered_data/", sample_name, "_cell_metadata.csv"))

cat("\n=== QC Pipeline Complete ===\n")
cat("Final cells:", ncol(seurat_obj), "\n")
cat("Final genes:", nrow(seurat_obj), "\n")
cat("Plots saved:", length(list.files("plots")), "files\n\n")

cat("Files created:\n")
cat(paste0("filtered_data/", sample_name, "_qc_filtered.rds - Clean Seurat object\n"))
cat(paste0("qc_metrics/", sample_name, "_QC_summary.csv - Comprehensive QC summary\n"))
cat(paste0("filtered_data/", sample_name, "_cell_metadata.csv - Cell-level metadata\n"))
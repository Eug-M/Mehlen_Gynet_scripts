
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
path <- "/home/eugenie-modolo/Documents/Gynet/scripts/Mehlen_Gynet_scripts/single_cell_QC_analysis"
setwd(path)
# setwd("~/GSE174609_scRNA/QC_analysis")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("qc_metrics", showWarnings = FALSE)
dir.create("filtered_data", showWarnings = FALSE)

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
sample_name <- "C1D1"
cellranger_output <- "/home/eugenie-modolo/Documents/Gynet/Samples/phase1_singlecell/"

# OPTION A: Load RAW matrix (recommanded)
counts <- Read10X(data.dir = file.path(cellranger_output, "C1D1_cellranger"))

# OPTION B: Load FILTERED matrix (alternative - skip Steps 4-5 if using this)
# counts <- Read10X(data.dir = file.path(cellranger_output, "filtered_feature_bc_matrix"))

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
seurat_obj$sample_id <- "C1D1"
seurat_obj$sra_id <- "SRRxxx"
seurat_obj$condition <- "C1D1"
seurat_obj$patient_id <- "Donor_C1D1"
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


#-----------------------------------------------
# STEP 4: Empty droplet detection - skipped for this dataset (we don't have the raw data before Cell Ranger)
#-----------------------------------------------

# Convert Seurat object to SingleCellExperiment for EmptyDrops
# We already loaded the raw matrix in Step 2
sce <- as.SingleCellExperiment(seurat_obj)

# Track initial droplet count for QC summary
initial_droplet_count <- ncol(sce)
initial_gene_count <- nrow(sce)

# Run EmptyDrops
set.seed(100)
empty_results <- emptyDrops(
  m = counts(sce),
  lower = 100,        # Droplets with <100 UMI assumed empty (for estimating ambient RNA)
  niters = 10000,     # Iterations for p-value calculation
  test.ambient = TRUE # Test if droplets differ from ambient RNA profile
)

# Identify cells (FDR < 0.01)
is_cell <- empty_results$FDR < 0.01
is_cell[is.na(is_cell)] <- FALSE  # Treat NA as empty (very low UMI droplets)
validated_barcodes <- colnames(sce)[is_cell]

# Brief feedback
cat("Droplets tested:", ncol(sce), "\n")
cat("Cells called:", sum(is_cell), 
    "(", round(sum(is_cell)/ncol(sce)*100, 1), "%)\n")
cat("Empty droplets removed:", sum(!is_cell), 
    "(", round(sum(!is_cell)/ncol(sce)*100, 1), "%)\n")

# Visual QC - the critical check!
empty_df <- data.frame(
  total_umi = colSums(counts(sce)),
  is_cell = is_cell
)

p1 <- ggplot(empty_df, aes(x = log10(total_umi + 1), fill = is_cell)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(
    values = c("TRUE" = "#2E86AB", "FALSE" = "#A23B72"),
    labels = c("Empty Droplet", "Cell"),
    name = "Classification"
  ) +
  labs(
    title = "EmptyDrops: Cell vs Empty Droplet Detection",
    subtitle = paste(sum(is_cell), "cells called from", ncol(sce), "total droplets"),
    x = "log10(UMI + 1)",
    y = "Number of Droplets"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/01_empty_droplets.png", p1, width = 10, height = 6, dpi = 300)

# NOW filter the Seurat object to keep only validated cells
# IMPORTANT: Save raw counts BEFORE filtering (needed for SoupX in Step 5)
raw_counts_all_droplets <- LayerData(seurat_obj, layer = "counts")

seurat_obj <- subset(seurat_obj, cells = validated_barcodes)

cat("After EmptyDrops:", ncol(seurat_obj), "cells retained\n")


#-----------------------------------------------
# STEP 5: Ambient RNA correction with SoupX - skipped for this dataset (we don't have the raw data before Cell Ranger)
#-----------------------------------------------

# Prepare SoupX channel
# tod (table of droplets) = raw counts from ALL droplets (saved in Step 4)
# toc (table of cells) = filtered counts from validated cells only
tod <- raw_counts_all_droplets
toc <- LayerData(seurat_obj, layer = "counts")
sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)

# Quick clustering for contamination estimation
temp_obj <- seurat_obj
temp_obj <- NormalizeData(temp_obj, verbose = FALSE)
temp_obj <- FindVariableFeatures(temp_obj, nfeatures = 2000, verbose = FALSE)
temp_obj <- ScaleData(temp_obj, verbose = FALSE)
temp_obj <- RunPCA(temp_obj, npcs = 30, verbose = FALSE)
temp_obj <- FindNeighbors(temp_obj, dims = 1:30, verbose = FALSE)
temp_obj <- FindClusters(temp_obj, resolution = 0.8, verbose = FALSE)
sc <- setClusters(sc, setNames(as.character(temp_obj$seurat_clusters), colnames(temp_obj)))

# Estimate contamination
sc <- tryCatch({
  autoEstCont(sc, verbose = FALSE)
}, error = function(e) {
  cat("Note: autoEstCont failed\n")
  sc$fit$rho <- NULL
  return(sc)
})

contamination_fraction <- sc$fit$rho

# Check contamination level
cat("Contamination fraction:", 
    ifelse(is.null(contamination_fraction), "NULL (very clean sample)", 
           paste0(round(contamination_fraction * 100, 2), "%")), "\n")

# Decision: Skip or correct based on contamination level
if (is.null(contamination_fraction) || contamination_fraction < 0.05) {
  cat("→ Contamination is NULL or <5% - skipping SoupX correction\n")
  cat("   Your sample is clean! Proceeding with original counts.\n")
} else {
  cat("→ Contamination is", round(contamination_fraction * 100, 2), "% (≥5%)\n")
  cat("   Applying SoupX correction...\n")
  
  # Apply correction
  suppressWarnings({
    corrected_counts <- adjustCounts(sc)
  })
  seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = corrected_counts)
  cat("   ✓ SoupX correction applied\n")
}



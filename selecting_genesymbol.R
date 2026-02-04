# Strategy for resolution of gene symbols (e.g. for GSEA) which correspond to several Ensembl id in counts_genes, in order of priority:
# 1. Prefer non-pseudogenes over pseudogenes
# 2. Prefer genes in keep_biotype over others
# 3. Prefer canonical chromosomes (1-22, X, Y) over patches/scaffolds
# 4. Prefer genes with RefSeq matches
# 5. If none of the above, sum counts for duplicate symbols



library(biomaRt)
library(dplyr)
library(stringr)



create_gene_symbol_mapping <- function(gene_info, keep_biotype, output_file) {
  
  cat("=== RESOLVING DUPLICATE GENE SYMBOLS ===\n\n")
  
  # Add resolution criteria
  gene_info_scored <- gene_info %>%
    mutate(
      is_not_pseudogene = !str_detect(gene_biotype, ".*pseudogene.*"),
      is_correct_biotype = gene_biotype %in% keep_biotype,
      is_main_chr = chromosome_name %in% c(1:22, "X", "Y", "MT"),
      has_refseq = !is.na(refseq_mrna) & refseq_mrna != "",
      n_duplicates = n(),
      .by = external_gene_name
    )
  
  # Separate unique and duplicated symbols
  unique_genes <- gene_info_scored %>%
    filter(n_duplicates == 1)
  
  dup_genes <- gene_info_scored %>%
    filter(n_duplicates > 1)
  
  cat("Unique gene symbols:", nrow(unique_genes), "\n")
  cat("Duplicated gene symbols:", length(unique(dup_genes$external_gene_name)), "\n")
  cat("Total genes with duplicated symbols:", nrow(dup_genes), "\n\n")
  
  # Process duplicates
  resolved_dups <- dup_genes %>%
    group_by(external_gene_name) %>%
    arrange(
      external_gene_name,
      desc(is_not_pseudogene),
      desc(is_correct_biotype),
      desc(has_refseq),
      desc(is_main_chr),
      ensembl_gene_id_version  # Tie-breaker: alphabetical
    ) %>%
    mutate(
      rank = row_number(),
      # Check if top-ranked gene is unique
      top_pseudogene = dplyr::first(is_not_pseudogene),
      # top_pseudogene = is_not_pseudogene,
      top_biotype = dplyr::first(is_correct_biotype),
      # top_biotype = is_correct_biotype,
      # top_biotype = Reduce(`|`, select(., is_correct_biotype)),
      top_refseq = dplyr::first(has_refseq),
      # top_refseq = has_refseq,
      top_chr = dplyr::first(is_main_chr),
      # top_chr = is_main_chr,
      # Count how many genes share the top score
      n_tied = sum(is_not_pseudogene == top_pseudogene & 
                     is_correct_biotype == top_biotype & 
                     has_refseq == top_refseq & 
                     is_main_chr == top_chr)
    ) %>%
    ungroup()
  
  # Split into resolved and still-tied
  uniquely_resolved <- resolved_dups %>%
    filter(n_tied == 1, rank == 1)
  
  still_tied <- resolved_dups %>%
    filter(n_tied > 1)
  
  cat("Resolved by criteria:", length(unique(uniquely_resolved$external_gene_name)), 
      "symbols\n")
  cat("Still tied (will be summed):", 
      length(unique(still_tied$external_gene_name)), "symbols\n\n")
  
  # Create mapping dataframe
  # For unique genes: one-to-one mapping
  mapping_unique <- unique_genes %>%
    mutate(
      final_gene_id = ensembl_gene_id_version,
      resolution_strategy = "unique",
      ensembl_ids_to_use = ensembl_gene_id_version,
      n_genes_merged = 1
    ) %>%
    dplyr::select(final_gene_id, external_gene_name, ensembl_ids_to_use, 
           resolution_strategy, n_genes_merged, gene_biotype, chromosome_name)
  
  # For uniquely resolved duplicates: one-to-one mapping
  mapping_resolved <- uniquely_resolved %>%
    filter(rank == 1) %>%
    mutate(
      final_gene_id = ensembl_gene_id_version,
      resolution_strategy = "criteria_resolved",
      ensembl_ids_to_use = ensembl_gene_id_version,
      n_genes_merged = 1
    ) %>%
    dplyr::select(final_gene_id, external_gene_name, ensembl_ids_to_use, 
           resolution_strategy, n_genes_merged, gene_biotype, chromosome_name)
  
  # For tied genes: create merged mapping
  mapping_merged <- still_tied %>%
    group_by(external_gene_name) %>%
    summarise(
      final_gene_id = paste0(dplyr::first(external_gene_name), "_merged"),
      ensembl_ids_to_use = paste(ensembl_gene_id_version, collapse = ";"),
      resolution_strategy = "summed",
      n_genes_merged = n(),
      gene_biotype = paste(unique(gene_biotype), collapse = ";"),
      chromosome_name = paste(unique(chromosome_name), collapse = ";"),
      .groups = "drop"
    )
  
  # Combine all mappings
  final_mapping <- bind_rows(
    mapping_unique,
    mapping_resolved,
    mapping_merged
  )
  
  # Print details of merged genes
  if (nrow(mapping_merged) > 0) {
    cat("Details of merged genes:\n")
    for (i in 1:nrow(mapping_merged)) {
      ids <- strsplit(mapping_merged$ensembl_ids_to_use[i], ";")[[1]]
      cat("  ", mapping_merged$external_gene_name[i], ": merging",
          mapping_merged$n_genes_merged[i], "genes (",
          paste(ids, collapse = ", "), ")\n")
    }
    cat("\n")
  }
  
  cat("=== FINAL SUMMARY ===\n")
  cat("Total gene symbols:", nrow(final_mapping), "\n")
  cat("Unique genes:", sum(final_mapping$resolution_strategy == "unique"), "\n")
  cat("Resolved duplicates:", sum(final_mapping$resolution_strategy == "criteria_resolved"), "\n")
  cat("Merged genes:", sum(final_mapping$resolution_strategy == "summed"), "\n\n")
  
  # Save mapping file
  write.csv(final_mapping, output_file, row.names = FALSE)
  cat("Mapping saved to:", output_file, "\n")
  
  return(final_mapping)
}



# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get all your Ensembl IDs
ensembl_ids <- rownames(counts_genes) # for 'ensembl_gene_id_version'
# ensembl_ids <- sub("\\..*", "", rownames(counts_genes)) # for 'ensembl_gene_id'

# Query biomaRt for transcript information
gene_info <- getBM(
  attributes = c(
    'ensembl_gene_id_version',        # Gene stable ID version
    # 'ensembl_gene_id',
    'external_gene_name',             # Gene name
    'gene_biotype',                   # Gene type
    # 'transcript_biotype',             # Transcript type
    'transcript_is_canonical',        # Ensembl Canonical
    'refseq_mrna',                    # RefSeq mRNA
    'chromosome_name'                 # Chromosome/scaffold name
    # 'transcript_mane_plus_clinical',  # RefSeq match transcript (MANE Plus Clinical) 
    # 'transcript_mane_select',         # RefSeq match transcript (MANE Select)
    # 'external_gene_source',           # Source of gene name
    # 'source',                         # Source (gene)
    # 'source_name',                    # Source name
    # 'name_1006',                      # GO term name
    # 'biogrid',                        # BioGRID Interaction data, The General Repository for Interaction Datasets ID
    # 'ccds',                           # CCDS ID
    # 'chembl',                         # ChEMBL ID
    # 'entrezgene_trans_name',          # EntrezGene transcript name ID
    # 'entrezgene_id'                  # NCBI gene (formerly Entrezgene) ID
    # 'uniprot_gn_id'                   # UniProtKB Gene Name ID
    # 'wikigene_name',                  # WikiGene name
    # 'smart'                           # SMART ID
  ),
  filters = 'ensembl_gene_id_version', #'ensembl_gene_id',
  values = ensembl_ids,
  mart = ensembl
)

# # Min & max number of transcripts per gene
# nber_duplicates <- list()
# for (id in ensembl_ids) { 
#   nber_duplicates <- cbind(nber_duplicates, nrow(gene_info[which(gene_info$ensembl_gene_id_version == id),])) 
# }
# min_value <- min(unlist(nber_duplicates)) # 1
# max_value <- max(unlist(nber_duplicates)) # 369


# We get duplicates, as we have one line per transcript, and not gene -> deduplication:
gene_info_summarized <- gene_info %>%
  group_by(ensembl_gene_id_version, external_gene_name, gene_biotype) %>%
  summarize(
    n_transcripts = n(),
    has_canonical = any(transcript_is_canonical == 1), # ils sont tous à TRUE
    refseq_mrna = paste(unique(refseq_mrna[!is.na(refseq_mrna)]), collapse = ", "),
    chromosome_name = paste(unique(chromosome_name[!is.na(chromosome_name)]), collapse = ", ")
  )
# gene_info_summarized$refseq_mrna <- ifelse(!is.na(gene_info_summarized$refseq_mrna) & gene_info_summarized$refseq_mrna != "", 1, 0)
gene_info_summarized <- data.frame(gene_info_summarized) # sinon, erreur "Can't supply `.by` when `.data` is a grouped data frame."
gene_info_summarized <- gene_info_summarized[which(gene_info_summarized$external_gene_name != ''), ]

keep_biotype <- c("protein_coding","IG_V_gene","IG_C_gene","IG_J_gene","IG_D_gene","TR_V_gene","TR_J_gene" ,"TR_C_gene","TR_D_gene")



mapping <- create_gene_symbol_mapping(
  gene_info_summarized, 
  keep_biotype, 
  output_file = "/home/eugenie-modolo/Documents/Reference_files/gene_symbol_mapping_Lapnet.csv"
)



# Function to apply the mapping to a count matrix
apply_gene_symbol_mapping <- function(count_matrix, mapping_file) {
  
  cat("=== APPLYING GENE SYMBOL MAPPING TO COUNT MATRIX ===\n\n")
  
  # Read mapping if it's a file path
  if (is.character(mapping_file)) {
    mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  } else {
    mapping <- mapping_file
  }
  
  cat("Mapping contains", nrow(mapping), "gene symbols\n")
  cat("Count matrix contains", nrow(count_matrix), "genes\n\n")
  
  # gene_name <- mapping$external_gene_name[i] 
  # Pre-compute which IDs exist in count_matrix 
  available_ids <- rownames(count_matrix)
  available_set <- as.list(setNames(rep(TRUE, length(available_ids)), available_ids))
  
  # Separate regular and merged
  is_merged <- grepl("_merged", mapping$final_gene_id)
  regular_genes <- mapping[!is_merged, ]
  merged_genes <- mapping[is_merged, ]
  
  # Process regular genes 
  regular_ids_to_keep <- regular_genes[match(available_ids, regular_genes$final_gene_id, nomatch = 0), 'final_gene_id']
  # regular_ids_to_keep <- regular_genes$final_gene_id[regular_genes$final_gene_id %in% available_ids]
  
  new_matrix <- count_matrix[regular_ids_to_keep, , drop = FALSE]
  rownames(new_matrix) <- regular_genes[match(available_ids, regular_genes$final_gene_id, nomatch = 0), 'external_gene_name']
  # rownames(new_matrix) <- regular_genes$external_gene_name[regular_genes$final_gene_id %in% available_ids]
  
  cat("Kept", length(regular_ids_to_keep), "regular genes\n")
  
  # Process merged genes (only if they exist)
  if (nrow(merged_genes) > 0) {
    cat("Processing", nrow(merged_genes), "merged genes...\n")
    
    # Pre-allocate matrix
    n_samples <- ncol(count_matrix)
    merged_matrix <- matrix(0, nrow = nrow(merged_genes), ncol = n_samples,
                            dimnames = list(merged_genes$external_gene_name,
                                            colnames(count_matrix)))
    
    # Loop only through merged genes 
    for (i in seq_len(nrow(merged_genes))) {
      original_ids <- unlist(strsplit(merged_genes$ensembl_ids_to_use[[i]], ';'))
      
      # Fast check which IDs exist
      existing_ids <- original_ids[original_ids %in% available_ids]
      
      if (length(existing_ids) > 0) {
        merged_matrix[i, ] <- colSums(count_matrix[existing_ids, , drop = FALSE])
      }
    }
    
    # Combine
    new_matrix <- rbind(new_matrix, merged_matrix)
  }
  
  cat("Output:", nrow(new_matrix), "genes x", ncol(new_matrix), "samples\n")
  
  return(new_matrix)
}


# Calling this function with the file name
count_matrix_resolved <- apply_gene_symbol_mapping(
  count_matrix, 
  mapping_file = "/home/eugenie-modolo/Documents/Reference_files/gene_symbol_mapping.csv"
)

# Calling this function with the mapping object in memory
count_matrix_resolved <- apply_gene_symbol_mapping(
  count_matrix, 
  mapping_file = mapping
)

# # Then create the dataframe for GSEA # obsolète, car j'ai déjà pris les gene symbols en rownames dans la fonction qui crée count_matrix_resolved
# count_for_gsea <- count_matrix_resolved
# # Get gene symbols from mapping
# gene_symbols <- mapping$external_gene_name[
#   match(rownames(count_for_gsea), mapping$final_gene_id)
# ]
# rownames(count_for_gsea) <- gene_symbols




# Testing the compatibility between the genes' biotypes from Ensembl and gtf annotation (from Gencode) 
gtf <- import('/home/eugenie-modolo/Documents/Reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gtf')
gtf_df <- as.data.frame(gtf)
for (ensbl_id in gene_info$gene_id) {
  if (gene_info[which(gene_info$gene_id == ensbl_id), 'gene_type'] != gtf_df[match(ensbl_id, gtf_df$gene_id, nomatch = 0), 'gene_type']) {
    print(c(ensbl_id, ' fonctionne pas : ', gene_info[which(gene_info$final_gene_id == ensbl_id), 'gene_type'], gtf_df[match(ensbl_id, gtf_df$gene_id, nomatch = 0), 'gene_type']))
  }
}
# nothing was printed out, so no incompatibility!



# # Filter for genes on chromosomes
# chromosomal_genes <- gene_info %>%
#   filter(grepl("^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$", chromosome_name)) %>% # to match chromosome names 1-22, X, Y, and MT
#   # filter(!grepl("^(GL|chrUn_|KI|random|scaffold)", chromosome_name, ignore.case = TRUE)) %>% # to exclude scaffolds, contigs and random chromosomes
#   distinct(ensembl_gene_id_version, external_gene_name, .keep_all = TRUE)
# 
# # Filter for genes with canonical transcripts
# canonical_genes <- gene_info %>%
#   filter(transcript_is_canonical == 1) %>%
#   distinct(ensembl_gene_id_version, external_gene_name, .keep_all = TRUE)
# 
# # Filter for genes with RefSeq matches
# refseq_genes <- gene_info %>%
#   filter(!is.na(refseq_mrna) & refseq_mrna != "") %>%
#   distinct(ensembl_gene_id_version, external_gene_name, .keep_all = TRUE)



# ## Use a pre-built annotation database
# 
# library(AnnotationHub)
# library(ensembldb)
# 
# # Access Ensembl database
# ah <- AnnotationHub()
# 
# # Query for human Ensembl database (adjust version as needed)
# query(ah, c("Homo sapiens", "EnsDb"))
# 
# # Get the database (use latest or specific version)
# edb <- ah[["AH113665"]]  # Example ID - check for latest
# 
# # Extract gene information
# genes_df <- genes(edb, return.type = "DataFrame")
# 
# # Get gene-to-symbol mapping with biotype
# gene_mapping <- as.data.frame(genes_df) %>%
#   select(gene_id, gene_name, gene_biotype, seq_name)
# 
# # Resolve duplicates
# unique_mapping <- gene_mapping %>%
#   group_by(gene_name) %>%
#   arrange(gene_name,
#           desc(gene_biotype == "protein_coding"),
#           gene_id) %>%
#   slice(1) %>%
#   ungroup()

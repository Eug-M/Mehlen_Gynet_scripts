<!--
SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# Mehlen_Gynet_scripts
R and Shell scripts from project GyNet, carried out at CRCL (Centre de Recherche en Cancérologie de Lyon) facilities, using data from biotech Netris pharma.


## Sub-project 1: bulk RNA-Seq data from clinical study (phase Ib) [work in progress]

The files cited hereafter were used to study the 168 samples obtained from 168 patients during the GyNet clinical trial. 

The goal of this trial was to assess if the addition of NP137 (Netrin-1 inhibitor) to the standard-of-care (chemo in arm B, immunotherapy in arm C, both in arm D) would increase survival in second- to sixth-line patients with either cervix or endometrial carcinomas. For this, samples where collected before receiving the treatment on 168 patients (67 with cervix carcinoma and 101 with endometrial carcinoma).

- selecting_genesymbol.R: script used in order to resolve the duplicates when one gene symbol (actual gene name) is linked to several ensembl gene IDs (ENSG*). I used biomart's web API in order to download, for the ensembl gene IDs present in my bulk RNA-seq results: their location on the genome, their biotype, their RefSeq match. For each conflict over a gene symbol, I then choose the best ensembl gene IDs using these information, and if no resolution is possible then I sum the ensembl gene IDs' expression.

- rapport_GyNet_1_QC.Rmd: QC (Quality Control) checks on the 168 samples (PCA, hierarchical clustering heatmaps, tumor purity), in order to see if there are samples behaving incoherently, and decide if they should be removed from further analyses or not. The secondary goal is to find which variables better explain the variable of interest of this study: survival.
As there is no consensus on the number of most-variable genes that should be used for the gene counts PCA, I decided to implement an elbow method on the variance (PC1 + PC2) explained by the number of top variable genes. 
As the ultimate goal will be to choose which variable(s) to include in the DESeq2 design (cf. files rapport_GyNet_3_choixDesign.Rmd & rapport_GyNet_4_DEGs.Rmd), I also performed PLS (Partial Least Square) and looked at the top variables contributing to the PCA variances, in order to have a first idea of which clinical variables most influence the gene expression in the samples.

- rapport_GyNet_2_survival.Rmd: study of the link between expression of Netrin-1, and its 15 main known receptors as well as global EMT (Epithelial-Mesenchymal Transition) levels, and the patients' survival (OS & PFS). The influence of the clinical trial's arm on the patient's survival is also studied at length (espacially arm A, which is the only one that did not receive NP137 treatment). This study notably uses Kaplan-Meier and ssGSEA statistical tests. 

- rapport_GyNet_3_choixDesign.Rmd: in addition to the PLS studies performed in the file rapport_GyNet_1_QC.Rmd, I performed here LRT (Likelihood Ratio Tests) with DESeq2 on the two different cancers (cervix and entrometrium), in order to automatically choose (i) the variable of interest as several clinical variables could be used to indicate survival; (ii) the best variables to include in the design to explain the noise in the data, so that more DEGs (Differentially-Expressed Genes) could be found.

- rapport_GyNet_4_DEGs.Rmd: study of the DEGs found for both cancers. I first choose the best DESeq2 design based on the heatmaps and Venn diagram, then perform cross-validation (re-run DESeq2 with the same design on all the samples minus one, for each sample) in order to keep only the DEGs validated by at least 75% of the samples. Finally I run the classical bio-informatics study on the DEGs found (ORA - Over-Enrichment Analysis, GSEA - Gene Set Enrichment Analysis).

- rapport_GyNet_UMIs.Rmd: comparison of the genes expressions with or without taking into account the UMIs (Unique Molecular Identifiers) inserted during the sequencing.

- using_puree.sh: Shell commands used to run the Python program PUREE to determine tumor purity in bulk RNA-Seq samples [PUREE: accurate pan-cancer tumor purity estimation from gene expression data](https://www.nature.com/articles/s42003-023-04764-8).

- commands_server.sh: Shell commands used on the lab's secure server (using slurm for distributed computing).


## Sub-project 2: single-cell RNA-Seq data from clinical study (phase Ia) [work in progress]

The files look into the single-cell RNA-Seq data from a patient with advanced endometrial carcinoma before and after treatment with NP137.
These data led to the publication of a paper in Nature, before the beginning of my work: [Netrin-1 blockade inhibits tumour growth and EMT features in endometrial cancer](https://www.nature.com/articles/s41586-023-06367-z).

- single_cell_*.R: study of the Netrin-1 expression, as well as it main receptors, on the single-cell RNA-Seq sequencing data. The different files are the various steps to process single-cell data. 


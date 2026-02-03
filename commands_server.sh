# Run nextflow without UMIs on Gnomic:
nextflow run -bg /home/modolo/nf-core-rnaseq_3.21.0/3_21_0/main.nf \
    --input /home/modolo/GyNet/Samples/samplesheet_GyNet.csv \
    --outdir /home/modolo/GyNet/Samples/results_test \
    --gtf /home/modolo/Reference_files/gencode.v49.primary_assembly.annotation.gtf \
    --fasta /home/modolo/Reference_files/GRCh38.primary_assembly.genome.fa \
    --aligner star_salmon \
    --extra_salmon_quant_args="--gcBias" \
    --min_mapped_reads 1 \
    --gencode \
    -profile server > log_nextflow.txt


# Run nextflow with UMIs on Gnomic:
nextflow run -bg /home/modolo/nf-core-rnaseq_3.21.0/3_21_0/main.nf \
    --input /home/modolo/GyNet_umis/230118_A00317_0594_BHN5TTDRX2/samplesheet_GyNet.csv \
    --outdir /home/modolo/GyNet_umis/230118_A00317_0594_BHN5TTDRX2/results_test \
    --gtf /home/modolo/Reference_files/gencode.v49.primary_assembly.annotation.gtf \
    --fasta /home/modolo/Reference_files/GRCh38.primary_assembly.genome.fa \
    --aligner star_salmon \
    --extra_salmon_quant_args="--gcBias" \
    --min_mapped_reads 1 \
    --gencode \
    --with_umi --umitools_extract_method "string" --umitools_bc_pattern "NNNN" --umitools_bc_pattern2 "NNNN" --save_umi_intermeds \
    -profile server > log_nextflow.txt


# Explore the bam files with IGV:
igv


# Look at the differences in one particular read found in IGV:
cd /Documents/Gynet/
samtools view ./results_rnaseq_NFcore/S05_002.markdup.sorted.bam | grep "A00317:594:HN5TTDRX2:1:2151:14895:9565"
samtools view ./results_rnaseq_NFcore_withUMIs/S05_002.umi_dedup.sorted.bam | grep "A00317:594:HN5TTDRX2:1:2151:14895:9565"


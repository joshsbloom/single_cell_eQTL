# A.vcf and  B.vcf pulled and extracted from
# /data/rrv2/genotyping/parents/cross_vcf/


# ..-2 are the samples with normal loading (#1 and #3)
vartrix --bam=/data/single_cell_eQTL/yeast2/3_ByxRM_51-_1-2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast2/3_ByxRM_51-_1-2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast2/3_ByxRM_51-_1-2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast2/3_ByxRM_51-_1-2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast2/3_ByxRM_51-_1-2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast2/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast2/A.vcf.recode.vcf \
    --threads=70 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info
#00:33:09 [INFO] Number of alignments evaluated: 3493207969
#00:33:09 [INFO] Number of alignments skipped due to low mapping quality: 2659971663
#00:33:09 [INFO] Number of alignments skipped due to not being primary: 0
#00:33:09 [INFO] Number of alignments skipped due to being duplicates: 0
#00:33:09 [INFO] Number of alignments skipped due to not being associated with a cell barcode: 32638005
#00:33:09 [INFO] Number of alignments skipped due to not being useful: 753988698
#00:33:09 [INFO] Number of alignments skipped due to not having a UMI: 5158
#00:33:09 [INFO] Number of VCF records skipped due to having invalid characters in the alternative haplotype: 0
#00:33:09 [INFO] Number of VCF records skipped due to being multi-allelic: 2667

vartrix --bam=/data/single_cell_eQTL/yeast2/1_2444_44_1-2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast2/1_2444_44_1-2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast2/1_2444_44_1-2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast2/1_2444_44_1-2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast2/1_2444_44_1-2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast2/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast2/B.vcf.recode.vcf \
    --threads=48 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

# ..- are the samples with overloading (#2 and #4)

vartrix --bam=/data/single_cell_eQTL/yeast2/2_2444_44-/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast2/2_2444_44-/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast2/2_2444_44-/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast2/2_2444_44-/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast2/2_2444_44-/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast2/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast2/B.vcf.recode.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=24 \
    --umi=TRUE \
    --log-level=info


vartrix --bam=/data/single_cell_eQTL/yeast2/4_ByxRM_51-/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast2/4_ByxRM_51-/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast2/4_ByxRM_51-/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast2/4_ByxRM_51-/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast2/4_ByxRM_51-/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast2/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast2/A.vcf.recode.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info



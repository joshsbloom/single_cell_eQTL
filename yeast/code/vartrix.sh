# A.vcf and  B.vcf pulled and extracted from
# /data/rrv2/genotyping/parents/cross_vcf/

#original BYxRM 480 set 
vartrix --bam=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_1/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_1/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_1/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_1/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_1/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/A.vcf.recode.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

vartrix --bam=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/00_BYxRM_480MatA_2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/A.vcf.recode.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info



# ..-2 are the samples with normal loading (#1 and #3)
vartrix --bam=/data/single_cell_eQTL/yeast/processed/03_ByxRM_51-_1-2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/03_ByxRM_51-_1-2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/03_ByxRM_51-_1-2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/03_ByxRM_51-_1-2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/03_ByxRM_51-_1-2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast2/reference/A.vcf.recode.vcf \
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

vartrix --bam=/data/single_cell_eQTL/yeast/processed/01_2444_44_1-2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/01_2444_44_1-2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/01_2444_44_1-2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/01_2444_44_1-2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/01_2444_44_1-2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/B.vcf.recode.vcf \
    --threads=48 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info



# ..- are the samples with overloading (#2 and #4)
vartrix --bam=/data/single_cell_eQTL/yeast/processed/02_2444_44-/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/02_2444_44-/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/02_2444_44-/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/02_2444_44-/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/02_2444_44-/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/processed/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/processed/reference/B.vcf.recode.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info


vartrix --bam=/data/single_cell_eQTL/yeast/processed/04_ByxRM_51-/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/04_ByxRM_51-/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/04_ByxRM_51-/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/04_ByxRM_51-/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/04_ByxRM_51-/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/processed/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/processed/A.vcf.recode.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

#generate merged vcf file for analysis of 4 parents 
#bcftools merge /data/single_cell_eQTL/yeast/reference/A.vcf.recode.vcf.gz /data/single_cell_eQTL/yeast/reference/B.vcf.recode.vcf.gz > BYxRMxYPS163xYJM145.vcf
#TO RUN

vartrix_linux --bam=/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/05_4-haploid-parental-strains/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/BYxRMxYPS163xYJM145.vcf \
    --threads=23 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

vartrix_linux --bam=/data/single_cell_eQTL/yeast/processed/07_2444-cross-1/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/07_2444-cross-1/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/07_2444-cross-1/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/07_2444-cross-1/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/07_2444-cross-1/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/B.vcf.recode.vcf \
    --threads=23 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

vartrix_linux --bam=/data/single_cell_eQTL/yeast/processed/07_2444-cross-2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/07_2444-cross-2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/07_2444-cross-2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/07_2444-cross-2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/07_2444-cross-2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/B.vcf.recode.vcf \
    --threads=23 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

vartrix --bam=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-1/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/parents.nostar.vcf \
    --threads=36 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

vartrix --bam=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/possorted_genome_bam.bam \
    --cell-barcodes=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/alt_counts.mtx \
    --ref-matrix=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/ref_counts.mtx \
    --out-variants=/data/single_cell_eQTL/yeast/processed/06_16-diploid-parental-strains-2/out_var.txt \
    --fasta=/data/single_cell_eQTL/yeast/reference/genome.fa \
    --vcf=/data/single_cell_eQTL/yeast/reference/parents.nostar.vcf \
    --threads=36 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info



./Local/vartrix/vartrix_linux --bam=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-1/possorted_genome_bam.bam \
    --cell-barcodes=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-1/filtered_feature_bc_matrix/barcodes.tsv \
       --out-matrix=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-1/alt_counts.mtx \
       --ref-matrix=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-1/ref_counts.mtx \
     --out-variants=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-1/out_var.txt \
            --fasta=/u/home/j/jsbloom/project-kruglyak/reference/genome.fa \
              --vcf=/u/home/j/jsbloom/project-kruglyak/reference/parents.nostar.vcf \
    --threads=24 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info

./Local/vartrix/vartrix_linux --bam=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-2/possorted_genome_bam.bam \
    --cell-barcodes=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-2/filtered_feature_bc_matrix/barcodes.tsv \
       --out-matrix=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-2/alt_counts.mtx \
       --ref-matrix=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-2/ref_counts.mtx \
     --out-variants=/u/home/j/jsbloom/project-kruglyak/processed/06_16-diploid-parental-strains-2/out_var.txt \
            --fasta=/u/home/j/jsbloom/project-kruglyak/reference/genome.fa \
              --vcf=/u/home/j/jsbloom/project-kruglyak/reference/parents.nostar.vcf \
    --threads=20 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info


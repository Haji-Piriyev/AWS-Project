# Variables
SRA = "SRR1972739"
REF_ID = "AF086833.2"
RESULTS_FOLDER = "results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR=f"{RESULTS_FOLDER}/snakemake"
BUCKET="haji-project1-bucket"
S3_PREFIX="ebola"
 
rule all:
    input: 
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{ALIGNED_DIR}/dedup.bam",
        f"{VARIANT_DIR}/raw_variants.vcf",
        f"{VARIANT_DIR}/filtered_variants.vcf",
        f"{ANNOTATED_DIR}/annotated_variants.vcf",
        f"{SNAKEMAKE_DIR}/.s3_upload_done"
       
 
rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        """
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """
 
rule download_reference:
    input:
        marker = rules.create_dirs.output.marker
    output:
        reference_fasta = f"{RAW_DIR}/reference.fasta"
    shell:
        """
        efetch -db nucleotide -id {REF_ID} -format fasta > {RAW_DIR}/reference.fasta
        samtools faidx {RAW_DIR}/reference.fasta
        bwa index {RAW_DIR}/reference.fasta
        gatk CreateSequenceDictionary -R {RAW_DIR}/reference.fasta -O {RAW_DIR}/reference.dict
        """
 
rule download_sra:
    input:
        marker = rules.create_dirs.output.marker
    output:
        sequence_sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        """
        prefetch {SRA} -O {RAW_DIR}
        """
 
rule extract_sequence:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra
    output:
        sequence_fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        fastqc -o {QC_DIR} {RAW_DIR}/{SRA}.fastq
        """

rule align_reads:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq
    output:
        sorted_dedupbam = f"{ALIGNED_DIR}/dedup.bam"

    shell:
        """
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {RAW_DIR}/reference.fasta {RAW_DIR}/{SRA}.fastq > {ALIGNED_DIR}/aligned.sam
        samtools view -b {ALIGNED_DIR}/aligned.sam > {ALIGNED_DIR}/aligned.bam  
        samtools sort {ALIGNED_DIR}/aligned.bam -o {ALIGNED_DIR}/aligned.sorted.bam
        gatk MarkDuplicates -I {ALIGNED_DIR}/aligned.sorted.bam -O {ALIGNED_DIR}/dedup.bam -M {ALIGNED_DIR}/dup_metrics.txt
        samtools index {ALIGNED_DIR}/dedup.bam
        """




rule call_variants:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        sorted_dedupbam = rules.align_reads.output.sorted_dedupbam
    output:
        raw_vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    shell:
        """
        gatk HaplotypeCaller -R {RAW_DIR}/reference.fasta -I {ALIGNED_DIR}/dedup.bam -O {VARIANT_DIR}/raw_variants.vcf
        """


rule filter_variants:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta,
        raw_vcf = rules.call_variants.output.raw_vcf
    output:
        filtered_vcf = f"{VARIANT_DIR}/filtered_variants.vcf"
    shell:
        """
        gatk VariantFiltration -R {RAW_DIR}/reference.fasta -V {VARIANT_DIR}/raw_variants.vcf -O {VARIANT_DIR}/filtered_variants.vcf --filter-expression "QD<2.0 || FS>60.0" --filter-name FILTER
        """


rule annotate_variants:
    input:
        filtered_vcf = rules.filter_variants.output.filtered_vcf
    output:
        annotated_vcf = f"{ANNOTATED_DIR}/annotated_variants.vcf"
    shell:
        """
        efetch -db nucleotide -id {REF_ID} -format genbank > {SNPEFF_DATA_DIR}/genes.gbk
        echo "reference_db.genome : reference_db" > {SNPEFF_DIR}/snpEff.config
        echo "reference_db.fa : {RAW_DIR}/reference.fasta" >> {SNPEFF_DIR}/snpEff.config
        echo "reference_db.genbank : {SNPEFF_DATA_DIR}/genes.gbk" >> {SNPEFF_DIR}/snpEff.config
        snpEff build -c {SNPEFF_DIR}/snpEff.config -genbank -v -noCheckProtein reference_db
        snpEff dump -c {SNPEFF_DIR}/snpEff.config reference_db > {SNPEFF_DIR}/snpEff_reference_db.txt 
        snpEff -c {SNPEFF_DIR}/snpEff.config -stats {SNPEFF_DIR}/snpEff.html reference_db {VARIANT_DIR}/filtered_variants.vcf > {ANNOTATED_DIR}/annotated_variants.vcf
        """



rule upload_s3:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta,
        sequence_sra = rules.download_sra.output.sequence_sra,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq,
        annotated_vcf = rules.annotate_variants.output.annotated_vcf
    output:
        marker = f"{SNAKEMAKE_DIR}/.s3_upload_done"
    run:
        import os
        import boto3
        s3 = boto3.client("s3")
 
        for root, dirs, files in os.walk(RESULTS_FOLDER):
            for file in files:
                local_file = os.path.join(root, file)
                relative_path = os.path.relpath(local_file, RESULTS_FOLDER)
                s3_key = os.path.join(S3_PREFIX, relative_path).replace("\\", "/")
 
                print(f"Uploading {local_file} to s3://{BUCKET}/{s3_key}")
                s3.upload_file(local_file, BUCKET, s3_key)
 
        with open(output.marker, "w") as f:
            f.write("Upload Complete!")
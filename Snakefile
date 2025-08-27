#!/usr/bin/env python3
"""
RNA-seq Analysis Pipeline using Snakemake
Processes paired-end RNA-seq data through QC, trimming, alignment, and quantification
"""

import os
import pandas as pd
from pathlib import Path

#############################################
#               CONFIGURATION               #
#############################################

# Reference genome files
REF_FASTA = "/path/to/maize_B73.fa"
REF_GTF = "/path/to/maize_B73.gtf"

# Pipeline parameters
READ_LENGTH = 150
SJDB_OVERHANG = READ_LENGTH - 1

# featureCounts parameters  
FC_STRANDED = 0        # 0=unstranded, 1=stranded, 2=reverse
FC_PAIRED = True       # Paired-end reads
FC_FEATURE = "exon"    # Feature type to count
FC_ATTR = "gene_id"    # Attribute for counting

# Resource allocation
THREADS = {
    "fastqc": 4,
    "fastp": 8, 
    "star_index": 16,
    "star_align": 12,
    "featurecounts": 8,
    "samtools": 4
}

# Output directories
DIRS = {
    "qc_raw": "results/qc/raw",
    "qc_trimmed": "results/qc/trimmed", 
    "trimmed": "results/trimmed",
    "star_index": "results/star_index",
    "aligned": "results/aligned",
    "counts": "results/counts",
    "multiqc": "results/multiqc"
}

# Sample definitions - 24 samples total
SAMPLES = {
    # Mock T1 samples
    "Mock_T1_R1": {"R1": "data/Mock_T1_R1_R1.fastq.gz", "R2": "data/Mock_T1_R1_R2.fastq.gz"},
    "Mock_T1_R2": {"R1": "data/Mock_T1_R2_R1.fastq.gz", "R2": "data/Mock_T1_R2_R2.fastq.gz"},
    "Mock_T1_R3": {"R1": "data/Mock_T1_R3_R1.fastq.gz", "R2": "data/Mock_T1_R3_R2.fastq.gz"},
    "Mock_T1_R4": {"R1": "data/Mock_T1_R4_R1.fastq.gz", "R2": "data/Mock_T1_R4_R2.fastq.gz"},
    # Mock T2 samples  
    "Mock_T2_R1": {"R1": "data/Mock_T2_R1_R1.fastq.gz", "R2": "data/Mock_T2_R1_R2.fastq.gz"},
    "Mock_T2_R2": {"R1": "data/Mock_T2_R2_R1.fastq.gz", "R2": "data/Mock_T2_R2_R2.fastq.gz"},
    "Mock_T2_R3": {"R1": "data/Mock_T2_R3_R1.fastq.gz", "R2": "data/Mock_T2_R3_R2.fastq.gz"},
    "Mock_T2_R4": {"R1": "data/Mock_T2_R4_R1.fastq.gz", "R2": "data/Mock_T2_R4_R2.fastq.gz"},
    # Treatment 1 T1 samples
    "Trt1_T1_R1": {"R1": "data/Trt1_T1_R1_R1.fastq.gz", "R2": "data/Trt1_T1_R1_R2.fastq.gz"},
    "Trt1_T1_R2": {"R1": "data/Trt1_T1_R2_R1.fastq.gz", "R2": "data/Trt1_T1_R2_R2.fastq.gz"},
    "Trt1_T1_R3": {"R1": "data/Trt1_T1_R3_R1.fastq.gz", "R2": "data/Trt1_T1_R3_R2.fastq.gz"},
    "Trt1_T1_R4": {"R1": "data/Trt1_T1_R4_R1.fastq.gz", "R2": "data/Trt1_T1_R4_R2.fastq.gz"},
    # Treatment 1 T2 samples
    "Trt1_T2_R1": {"R1": "data/Trt1_T2_R1_R1.fastq.gz", "R2": "data/Trt1_T2_R1_R2.fastq.gz"},
    "Trt1_T2_R2": {"R1": "data/Trt1_T2_R2_R1.fastq.gz", "R2": "data/Trt1_T2_R2_R2.fastq.gz"},
    "Trt1_T2_R3": {"R1": "data/Trt1_T2_R3_R1.fastq.gz", "R2": "data/Trt1_T2_R3_R2.fastq.gz"},
    "Trt1_T2_R4": {"R1": "data/Trt1_T2_R4_R1.fastq.gz", "R2": "data/Trt1_T2_R4_R2.fastq.gz"},
    # Treatment 2 T1 samples
    "Trt2_T1_R1": {"R1": "data/Trt2_T1_R1_R1.fastq.gz", "R2": "data/Trt2_T1_R1_R2.fastq.gz"},
    "Trt2_T1_R2": {"R1": "data/Trt2_T1_R2_R1.fastq.gz", "R2": "data/Trt2_T1_R2_R2.fastq.gz"},
    "Trt2_T1_R3": {"R1": "data/Trt2_T1_R3_R1.fastq.gz", "R2": "data/Trt2_T1_R3_R2.fastq.gz"},
    "Trt2_T1_R4": {"R1": "data/Trt2_T1_R4_R1.fastq.gz", "R2": "data/Trt2_T1_R4_R2.fastq.gz"},
    # Treatment 2 T2 samples
    "Trt2_T2_R1": {"R1": "data/Trt2_T2_R1_R1.fastq.gz", "R2": "data/Trt2_T2_R1_R2.fastq.gz"},
    "Trt2_T2_R2": {"R1": "data/Trt2_T2_R2_R1.fastq.gz", "R2": "data/Trt2_T2_R2_R2.fastq.gz"},
    "Trt2_T2_R3": {"R1": "data/Trt2_T2_R3_R1.fastq.gz", "R2": "data/Trt2_T2_R3_R2.fastq.gz"},
    "Trt2_T2_R4": {"R1": "data/Trt2_T2_R4_R1.fastq.gz", "R2": "data/Trt2_T2_R4_R2.fastq.gz"}
}

#############################################
#            PIPELINE SETUP                #
#############################################

# Extract sample names
SAMPLE_NAMES = list(SAMPLES.keys())

# Create output directories
for directory in DIRS.values():
    Path(directory).mkdir(parents=True, exist_ok=True)

#############################################
#               RULES                       #
#############################################

rule all:
    """Target rule - defines all final outputs"""
    input:
        # STAR index
        f"{DIRS['star_index']}/SAindex",
        # Quality control reports
        expand(f"{DIRS['qc_raw']}/{{sample}}_R1_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{DIRS['qc_raw']}/{{sample}}_R2_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{DIRS['qc_trimmed']}/{{sample}}_R1_trimmed_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{DIRS['qc_trimmed']}/{{sample}}_R2_trimmed_fastqc.html", sample=SAMPLE_NAMES),
        # Trimmed reads
        expand(f"{DIRS['trimmed']}/{{sample}}_R1_trimmed.fastq.gz", sample=SAMPLE_NAMES),
        expand(f"{DIRS['trimmed']}/{{sample}}_R2_trimmed.fastq.gz", sample=SAMPLE_NAMES),
        # Aligned reads
        expand(f"{DIRS['aligned']}/{{sample}}.bam", sample=SAMPLE_NAMES),
        expand(f"{DIRS['aligned']}/{{sample}}.bam.bai", sample=SAMPLE_NAMES),
        # Gene counts
        f"{DIRS['counts']}/gene_counts.txt",
        f"{DIRS['counts']}/count_matrix.tsv",
        # MultiQC report
        f"{DIRS['multiqc']}/multiqc_report.html"

rule build_star_index:
    """Build STAR genome index"""
    input:
        fasta=REF_FASTA,
        gtf=REF_GTF
    output:
        f"{DIRS['star_index']}/SAindex"
    params:
        index_dir=DIRS['star_index']
    threads: THREADS['star_index']
    resources:
        mem_mb=32000
    conda: "envs/star.yaml"
    log: "logs/star_index.log"
    shell:
        """
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {params.index_dir} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {SJDB_OVERHANG} \
             --limitGenomeGenerateRAM 30000000000 2> {log}
        """

rule fastqc_raw:
    """Quality control of raw reads"""
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        html_r1=f"{DIRS['qc_raw']}/{{sample}}_R1_fastqc.html",
        html_r2=f"{DIRS['qc_raw']}/{{sample}}_R2_fastqc.html",
        zip_r1=f"{DIRS['qc_raw']}/{{sample}}_R1_fastqc.zip",
        zip_r2=f"{DIRS['qc_raw']}/{{sample}}_R2_fastqc.zip"
    params:
        outdir=DIRS['qc_raw']
    threads: THREADS['fastqc']
    conda: "envs/qc.yaml"
    log: "logs/fastqc_raw_{sample}.log"
    shell:
        """
        # Create properly named symlinks for FastQC
        ln -sf $(readlink -f {input.r1}) {params.outdir}/{wildcards.sample}_R1.fastq.gz
        ln -sf $(readlink -f {input.r2}) {params.outdir}/{wildcards.sample}_R2.fastq.gz
        
        fastqc --threads {threads} --outdir {params.outdir} \
               {params.outdir}/{wildcards.sample}_R1.fastq.gz \
               {params.outdir}/{wildcards.sample}_R2.fastq.gz 2> {log}
               
        # Clean up symlinks
        rm {params.outdir}/{wildcards.sample}_R1.fastq.gz
        rm {params.outdir}/{wildcards.sample}_R2.fastq.gz
        """

rule trim_reads:
    """Trim adapters and low-quality bases with fastp"""
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        r1=f"{DIRS['trimmed']}/{{sample}}_R1_trimmed.fastq.gz",
        r2=f"{DIRS['trimmed']}/{{sample}}_R2_trimmed.fastq.gz",
        html=f"{DIRS['trimmed']}/{{sample}}_fastp.html",
        json=f"{DIRS['trimmed']}/{{sample}}_fastp.json"
    threads: THREADS['fastp']
    conda: "envs/qc.yaml"
    log: "logs/fastp_{sample}.log"
    shell:
        """
        fastp --thread {threads} \
              --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1} --out2 {output.r2} \
              --detect_adapter_for_pe \
              --qualified_quality_phred 20 \
              --unqualified_percent_limit 20 \
              --length_required 30 \
              --cut_front \
              --cut_tail \
              --cut_mean_quality 20 \
              --html {output.html} \
              --json {output.json} 2> {log}
        """

rule fastqc_trimmed:
    """Quality control of trimmed reads"""
    input:
        r1=f"{DIRS['trimmed']}/{{sample}}_R1_trimmed.fastq.gz",
        r2=f"{DIRS['trimmed']}/{{sample}}_R2_trimmed.fastq.gz"
    output:
        html_r1=f"{DIRS['qc_trimmed']}/{{sample}}_R1_trimmed_fastqc.html",
        html_r2=f"{DIRS['qc_trimmed']}/{{sample}}_R2_trimmed_fastqc.html",
        zip_r1=f"{DIRS['qc_trimmed']}/{{sample}}_R1_trimmed_fastqc.zip",
        zip_r2=f"{DIRS['qc_trimmed']}/{{sample}}_R2_trimmed_fastqc.zip"
    params:
        outdir=DIRS['qc_trimmed']
    threads: THREADS['fastqc']
    conda: "envs/qc.yaml"
    log: "logs/fastqc_trimmed_{sample}.log"
    shell:
        """
        fastqc --threads {threads} --outdir {params.outdir} \
               {input.r1} {input.r2} 2> {log}
        """

rule star_align:
    """Align reads with STAR"""
    input:
        r1=f"{DIRS['trimmed']}/{{sample}}_R1_trimmed.fastq.gz",
        r2=f"{DIRS['trimmed']}/{{sample}}_R2_trimmed.fastq.gz",
        index=f"{DIRS['star_index']}/SAindex"
    output:
        bam=f"{DIRS['aligned']}/{{sample}}.bam",
        log_final=f"{DIRS['aligned']}/{{sample}}_Log.final.out",
        log_out=f"{DIRS['aligned']}/{{sample}}_Log.out",
        log_progress=f"{DIRS['aligned']}/{{sample}}_Log.progress.out"
    params:
        index_dir=DIRS['star_index'],
        out_prefix=f"{DIRS['aligned']}/{{sample}}_",
        tmp_bam=f"{DIRS['aligned']}/{{sample}}_Aligned.sortedByCoord.out.bam"
    threads: THREADS['star_align']
    resources:
        mem_mb=32000
    conda: "envs/star.yaml"
    log: "logs/star_align_{sample}.log"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {params.index_dir} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMstrandField intronMotif \
             --outFilterIntronMotifs RemoveNoncanonical \
             --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
             --quantMode GeneCounts \
             --twopassMode Basic \
             --outFilterMultimapNmax 20 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             --outFilterMismatchNmax 999 \
             --outFilterMismatchNoverReadLmax 0.04 \
             --alignIntronMin 20 \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000 2> {log}
        
        # Move the output BAM file
        mv {params.tmp_bam} {output.bam}
        """

rule index_bam:
    """Index BAM files"""
    input:
        f"{DIRS['aligned']}/{{sample}}.bam"
    output:
        f"{DIRS['aligned']}/{{sample}}.bam.bai"
    threads: THREADS['samtools']
    conda: "envs/samtools.yaml"
    log: "logs/index_bam_{sample}.log"
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """

rule count_features:
    """Count reads mapping to genes with featureCounts"""
    input:
        bams=expand(f"{DIRS['aligned']}/{{sample}}.bam", sample=SAMPLE_NAMES),
        gtf=REF_GTF
    output:
        counts=f"{DIRS['counts']}/gene_counts.txt",
        summary=f"{DIRS['counts']}/gene_counts.txt.summary"
    params:
        strand_param=f"-s {FC_STRANDED}",
        paired_param="-p" if FC_PAIRED else ""
    threads: THREADS['featurecounts']
    conda: "envs/subread.yaml"
    log: "logs/featurecounts.log"
    shell:
        """
        featureCounts -T {threads} \
                      -a {input.gtf} \
                      -o {output.counts} \
                      -t {FC_FEATURE} \
                      -g {FC_ATTR} \
                      {params.strand_param} \
                      {params.paired_param} \
                      -B -C \
                      {input.bams} 2> {log}
        """

rule create_count_matrix:
    """Create a clean count matrix from featureCounts output"""
    input:
        f"{DIRS['counts']}/gene_counts.txt"
    output:
        f"{DIRS['counts']}/count_matrix.tsv"
    run:
        # Read featureCounts output
        df = pd.read_csv(input[0], sep='\t', comment='#', low_memory=False)
        
        # Extract gene info and count columns
        gene_info = df[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']]
        count_cols = [col for col in df.columns if col.endswith('.bam')]
        
        # Create count matrix with clean sample names
        count_data = df[['Geneid'] + count_cols].copy()
        
        # Clean up column names (remove path and .bam extension)
        new_columns = ['gene_id']
        for col in count_cols:
            sample_name = os.path.basename(col).replace('.bam', '')
            new_columns.append(sample_name)
        
        count_data.columns = new_columns
        
        # Save count matrix
        count_data.to_csv(output[0], sep='\t', index=False)
        
        # Also save gene info
        gene_info.to_csv(f"{DIRS['counts']}/gene_info.tsv", sep='\t', index=False)

rule multiqc:
    """Generate MultiQC report"""
    input:
        # Raw FastQC
        expand(f"{DIRS['qc_raw']}/{{sample}}_R1_fastqc.zip", sample=SAMPLE_NAMES),
        expand(f"{DIRS['qc_raw']}/{{sample}}_R2_fastqc.zip", sample=SAMPLE_NAMES),
        # Trimmed FastQC  
        expand(f"{DIRS['qc_trimmed']}/{{sample}}_R1_trimmed_fastqc.zip", sample=SAMPLE_NAMES),
        expand(f"{DIRS['qc_trimmed']}/{{sample}}_R2_trimmed_fastqc.zip", sample=SAMPLE_NAMES),
        # fastp reports
        expand(f"{DIRS['trimmed']}/{{sample}}_fastp.json", sample=SAMPLE_NAMES),
        # STAR logs
        expand(f"{DIRS['aligned']}/{{sample}}_Log.final.out", sample=SAMPLE_NAMES),
        # featureCounts summary
        f"{DIRS['counts']}/gene_counts.txt.summary"
    output:
        f"{DIRS['multiqc']}/multiqc_report.html"
    params:
        outdir=DIRS['multiqc'],
        search_dirs=' '.join(DIRS.values())
    conda: "envs/multiqc.yaml"
    log: "logs/multiqc.log"
    shell:
        """
        multiqc {params.search_dirs} \
                --outdir {params.outdir} \
                --force \
                --verbose \
                --config multiqc_config.yaml 2> {log}
        """

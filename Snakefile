#############################################
#               USER SETTINGS               #
#############################################

# ---- Reference (Maize B73) ----
REF_FASTA  = "/path/to/maize_B73.fa"
REF_GTF    = "/path/to/maize_B73.gtf"
STAR_INDEX = "ref/star_index"   # will be created if missing
READ_LENGTH = 150               # typical 150bp; sjdbOverhang = 149

# ---- Threads per step ----
THREADS = {
    "fastqc": 2,
    "fastp": 4,
    "star_index": 16,
    "star_align": 8,
    "samtools": 4,
    "featurecounts": 8,
}

# ---- featureCounts parameters ----
FC_STRANDED = 0        # 0: unstranded, 1: stranded, 2: reverse
FC_PAIRED   = True     # Paired-end data?
FC_TYPE     = "exon"   # -t exon
FC_ID_ATTR  = "gene_id"  # -g gene_id

# ---- Output directories ----
OUT_QC_RAW  = "results/fastqc/raw"
OUT_QC_TRIM = "results/fastqc/trimmed"
OUT_FASTP   = "results/fastp"
OUT_STAR    = "results/star"
OUT_COUNTS  = "results/counts"
OUT_MQC     = "results/multiqc"

# ---- Samples (24 total; replace paths with your own) ----
# Use any filenames; only sample IDs must be unique.
SAMPLES = {
    # Mock T1
    "Mock_T1_R1": {"r1": "data/Mock_T1_R1_R1.fastq.gz", "r2": "data/Mock_T1_R1_R2.fastq.gz"},
    "Mock_T1_R2": {"r1": "data/Mock_T1_R2_R1.fastq.gz", "r2": "data/Mock_T1_R2_R2.fastq.gz"},
    "Mock_T1_R3": {"r1": "data/Mock_T1_R3_R1.fastq.gz", "r2": "data/Mock_T1_R3_R2.fastq.gz"},
    "Mock_T1_R4": {"r1": "data/Mock_T1_R4_R1.fastq.gz", "r2": "data/Mock_T1_R4_R2.fastq.gz"},
    # Mock T2
    "Mock_T2_R1": {"r1": "data/Mock_T2_R1_R1.fastq.gz", "r2": "data/Mock_T2_R1_R2.fastq.gz"},
    "Mock_T2_R2": {"r1": "data/Mock_T2_R2_R1.fastq.gz", "r2": "data/Mock_T2_R2_R2.fastq.gz"},
    "Mock_T2_R3": {"r1": "data/Mock_T2_R3_R1.fastq.gz", "r2": "data/Mock_T2_R3_R2.fastq.gz"},
    "Mock_T2_R4": {"r1": "data/Mock_T2_R4_R1.fastq.gz", "r2": "data/Mock_T2_R4_R2.fastq.gz"},
    # treatment-1 T1
    "Trt1_T1_R1": {"r1": "data/Trt1_T1_R1_R1.fastq.gz", "r2": "data/Trt1_T1_R1_R2.fastq.gz"},
    "Trt1_T1_R2": {"r1": "data/Trt1_T1_R2_R1.fastq.gz", "r2": "data/Trt1_T1_R2_R2.fastq.gz"},
    "Trt1_T1_R3": {"r1": "data/Trt1_T1_R3_R1.fastq.gz", "r2": "data/Trt1_T1_R3_R2.fastq.gz"},
    "Trt1_T1_R4": {"r1": "data/Trt1_T1_R4_R1.fastq.gz", "r2": "data/Trt1_T1_R4_R2.fastq.gz"},
    # treatment-1 T2
    "Trt1_T2_R1": {"r1": "data/Trt1_T2_R1_R1.fastq.gz", "r2": "data/Trt1_T2_R1_R2.fastq.gz"},
    "Trt1_T2_R2": {"r1": "data/Trt1_T2_R2_R1.fastq.gz", "r2": "data/Trt1_T2_R2_R2.fastq.gz"},
    "Trt1_T2_R3": {"r1": "data/Trt1_T2_R3_R1.fastq.gz", "r2": "data/Trt1_T2_R3_R2.fastq.gz"},
    "Trt1_T2_R4": {"r1": "data/Trt1_T2_R4_R1.fastq.gz", "r2": "data/Trt1_T2_R4_R2.fastq.gz"},
    # treatment-2 T1
    "Trt2_T1_R1": {"r1": "data/Trt2_T1_R1_R1.fastq.gz", "r2": "data/Trt2_T1_R1_R2.fastq.gz"},
    "Trt2_T1_R2": {"r1": "data/Trt2_T1_R2_R1.fastq.gz", "r2": "data/Trt2_T1_R2_R2.fastq.gz"},
    "Trt2_T1_R3": {"r1": "data/Trt2_T1_R3_R1.fastq.gz", "r2": "data/Trt2_T1_R3_R2.fastq.gz"},
    "Trt2_T1_R4": {"r1": "data/Trt2_T1_R4_R1.fastq.gz", "r2": "data/Trt2_T1_R4_R2.fastq.gz"},
    # treatment-2 T2
    "Trt2_T2_R1": {"r1": "data/Trt2_T2_R1_R1.fastq.gz", "r2": "data/Trt2_T2_R1_R2.fastq.gz"},
    "Trt2_T2_R2": {"r1": "data/Trt2_T2_R2_R1.fastq.gz", "r2": "data/Trt2_T2_R2_R2.fastq.gz"},
    "Trt2_T2_R3": {"r1": "data/Trt2_T2_R3_R1.fastq.gz", "r2": "data/Trt2_T2_R3_R2.fastq.gz"},
    "Trt2_T2_R4": {"r1": "data/Trt2_T2_R4_R1.fastq.gz", "r2": "data/Trt2_T2_R4_R2.fastq.gz"},
}

#############################################
#            PIPELINE (NO CONFIG)           #
#############################################

import os
from pathlib import Path

SAMPLE_IDS = sorted(SAMPLES.keys())
SJDB_OVERH = max(1, READ_LENGTH - 1)

# ensure directories exist
for d in [OUT_QC_RAW, OUT_QC_TRIM, OUT_FASTP, OUT_STAR, OUT_COUNTS, OUT_MQC, STAR_INDEX]:
    Path(d).mkdir(parents=True, exist_ok=True)

def trimmed_r1(wc): return f"{OUT_FASTP}/{wc.sample}_R1.trimmed.fastq.gz"
def trimmed_r2(wc): return f"{OUT_FASTP}/{wc.sample}_R2.trimmed.fastq.gz"
def bam_path(wc):   return f"{OUT_STAR}/{wc.sample}.bam"
def bai_path(wc):   return f"{OUT_STAR}/{wc.sample}.bam.bai"

rule all:
    input:
        f"{STAR_INDEX}/SAindex",
        expand(trimmed_r1, sample=SAMPLE_IDS),
        expand(trimmed_r2, sample=SAMPLE_IDS),
        expand(bam_path,   sample=SAMPLE_IDS),
        expand(bai_path,   sample=SAMPLE_IDS),
        expand(f"{OUT_QC_RAW}" + "/{sample}.fastqc.done", sample=SAMPLE_IDS),
        expand(f"{OUT_QC_TRIM}"+ "/{sample}.fastqc.done", sample=SAMPLE_IDS),
        f"{OUT_COUNTS}/featurecounts.txt",
        f"{OUT_COUNTS}/counts_matrix.tsv",
        f"{OUT_MQC}/multiqc_report.html"

# STAR genome index
rule star_index:
    output:
        f"{STAR_INDEX}/SAindex"
    input:
        fasta=REF_FASTA,
        gtf=REF_GTF
    threads: THREADS["star_index"]
    conda: "envs/align.yaml"
    shell:
        r"""
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {STAR_INDEX} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {SJDB_OVERH}
        test -s {output} || (echo "STAR index missing {output}" && exit 1)
        """

# Raw FastQC (marker files)
rule fastqc_raw:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["r1"],
        r2=lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        touch(f"{OUT_QC_RAW}" + "/{sample}.fastqc.done")
    threads: THREADS["fastqc"]
    conda: "envs/qc.yaml"
    shell:
        r"""
        fastqc -t {threads} -o {OUT_QC_RAW} {input.r1} {input.r2}
        """

# Trimming
rule fastp:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["r1"],
        r2=lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        r1=formatted("{OUT_FASTP}" + "/{sample}_R1.trimmed.fastq.gz"),
        r2=formatted("{OUT_FASTP}" + "/{sample}_R2.trimmed.fastq.gz"),
        html=formatted("{OUT_FASTP}" + "/{sample}.fastp.html"),
        json=formatted("{OUT_FASTP}" + "/{sample}.fastp.json")
    threads: THREADS["fastp"]
    conda: "envs/qc.yaml"
    shell:
        r"""
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --thread {threads} \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --length_required 30 \
          --html {output.html} \
          --json {output.json}
        """

# FastQC on trimmed reads (marker)
rule fastqc_trimmed:
    input:
        r1=trimmed_r1,
        r2=trimmed_r2
    output:
        touch(f"{OUT_QC_TRIM}" + "/{sample}.fastqc.done")
    threads: THREADS["fastqc"]
    conda: "envs/qc.yaml"
    shell:
        r"""
        fastqc -t {threads} -o {OUT_QC_TRIM} {input.r1} {input.r2}
        """

# STAR alignment → sorted BAM + index
rule star_align:
    input:
        idx=f"{STAR_INDEX}/SAindex",
        r1=trimmed_r1,
        r2=trimmed_r2
    output:
        bam=bam_path,
        bai=bai_path
    threads: THREADS["star_align"]
    conda: "envs/align.yaml"
    params:
        prefix=lambda wc: f"{OUT_STAR}/{wc.sample}."
    shell:
        r"""
        STAR --runThreadN {threads} \
             --genomeDir {STAR_INDEX} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --twopassMode Basic \
             --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample}
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# featureCounts (one run over all BAMs)
rule featurecounts:
    input:
        bams=expand(bam_path, sample=SAMPLE_IDS),
        gtf=REF_GTF
    output:
        txt=f"{OUT_COUNTS}/featurecounts.txt",
        summary=f"{OUT_COUNTS}/featurecounts.summary"
    threads: THREADS["featurecounts"]
    conda: "envs/align.yaml"
    shell:
        r"""
        featureCounts \
          -T {threads} \
          -a {input.gtf} \
          -o {output.txt} \
          -t {FC_TYPE} \
          -g {FC_ID_ATTR} \
          {'-p -B -C' if FC_PAIRED else ''} \
          {'-s ' + str(FC_STRANDED) if FC_STRANDED in [0,1,2] else ''} \
          {input.bams}
        mv {output.txt}.summary {output.summary}
        """

# Make a tidy counts matrix (gene_id + sample columns)
rule make_counts_matrix:
    input:
        f"{OUT_COUNTS}/featurecounts.txt"
    output:
        f"{OUT_COUNTS}/counts_matrix.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t", comment="#")
        keep_cols = ["Geneid"] + [c for c in df.columns if c.endswith(".bam")]
        df = df[keep_cols]
        bam_to_sample = {f"{OUT_STAR}/{s}.bam": s for s in SAMPLE_IDS}
        new_cols = ["gene_id"] + [bam_to_sample.get(c, Path(c).stem.replace(".bam","")) for c in df.columns[1:]]
        df.columns = new_cols
        df.to_csv(output[0], sep="\t", index=False)

# MultiQC summary
rule multiqc:
    input:
        expand(f"{OUT_QC_RAW}" + "/{sample}.fastqc.done",  sample=SAMPLE_IDS),
        expand(f"{OUT_QC_TRIM}"+ "/{sample}.fastqc.done",  sample=SAMPLE_IDS),
        expand(trimmed_r1, sample=SAMPLE_IDS),
        expand(trimmed_r2, sample=SAMPLE_IDS),
        expand(bam_path,   sample=SAMPLE_IDS),
        f"{OUT_COUNTS}/featurecounts.txt"
    output:
        html=f"{OUT_MQC}/multiqc_report.html"
    threads: 2
    conda: "envs/qc.yaml"
    shell:
        r"""
        multiqc -o {OUT_MQC} .
        """

# 🧬 Snakemake tutorial for a simple RNA-seq workflow

Welcome! 👋
This repo contains a **reproducible RNA-seq pipeline** built with [Snakemake](https://snakemake.readthedocs.io/). It takes you from **raw FASTQ files → quality control → trimming → alignment → read counts → tidy matrix** — all in one place.

The example here uses **Maize B73** as the reference genome and 24 RNA-seq samples (mock vs. two treatments, two timepoints, four biological replicates each). But you can easily swap in your own samples or organism.

<img src="flowchart.png" alt="Workflow" width="60%"/>

## 🌟 What This Pipeline Does

RNA sequencing (RNA-seq) is like taking a snapshot of all the genes being actively used in your cells at a specific moment. Imagine your genome as a massive library, and RNA-seq tells us which "books" (genes) are being "read" (expressed) and how frequently.

This pipeline takes your raw RNA-seq data and transforms it into meaningful gene expression counts through several carefully orchestrated steps—think of it as a sophisticated assembly line for genomic data.

## 🎯 What You'll Learn

By the end of this tutorial, you'll understand how to:
- Set up a reproducible RNA-seq analysis workflow
- Quality control your sequencing data
- Align reads to a reference genome
- Count gene expression levels
- Generate comprehensive analysis reports

## 📋 Prerequisites

**What you need to know:**
- Basic command line navigation (`cd`, `ls`, `mkdir`)
- What RNA-seq data looks like (FASTQ files)
- Basic understanding of bioinformatics concepts

**What you need installed:**
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/) (for managing software)
- [Snakemake](https://snakemake.readthedocs.io/) (our workflow manager)

```bash
# Install Snakemake
conda install -c conda-forge -c bioconda snakemake
```

## 📁 Setting Up Your Project

First, let's create a well-organized project structure. Think of this as setting up your laboratory workspace:

```bash
# Create your main project directory
mkdir rnaseq_analysis
cd rnaseq_analysis

# Create subdirectories
mkdir -p data envs logs
```

Your project will look like this:
```
rnaseq_analysis/
├── Snakefile                 # The main pipeline (your recipe book)
├── multiqc_config.yaml      # Report customization settings
├── envs/                    # Software environment definitions
├── data/                    # Your raw FASTQ files go here
├── logs/                    # Will store detailed run logs
└── results/                 # Will contain all analysis outputs
```

## 🔧 Understanding the Pipeline Components

### The Snakefile: Your Analysis Recipe

The `Snakefile` is like a detailed cookbook that tells the computer exactly how to process your data. Each "rule" is a recipe step:

1. **Quality Control (FastQC)**: Examines your raw data quality—like inspecting ingredients before cooking
2. **Read Trimming (fastp)**: Removes low-quality bases and adapters—like cleaning and prepping ingredients
3. **Genome Indexing (STAR)**: Prepares the reference genome for fast searching—like organizing your spice rack
4. **Read Alignment (STAR)**: Maps your reads to the genome—like following GPS directions
5. **Gene Counting (featureCounts)**: Tallies reads per gene—like counting votes in an election
6. **Report Generation (MultiQC)**: Creates a beautiful summary report—like writing up your experimental results

### Environment Files: Software Management Made Easy

Instead of manually installing dozens of bioinformatics tools (and dealing with version conflicts), we use conda environments. Think of each environment file as a shopping list for specific software packages.

## 🚀 Step-by-Step Setup

### Step 1: Download the Pipeline Files

Save the main pipeline as `Snakefile` in your project directory (copy from the artifact above).

### Step 2: Create Environment Files

Create the `envs/` directory and save each environment file:

**envs/star.yaml** - For genome indexing and read alignment:
```yaml
name: star
channels:
  - conda-forge
  - bioconda
dependencies:
  - star=2.7.10a
```

**envs/qc.yaml** - For quality control and read trimming:
```yaml
name: qc
channels:
  - conda-forge
  - bioconda
dependencies:
  - fastqc=0.11.9
  - fastp=0.23.2
```

**envs/samtools.yaml** - For BAM file manipulation:
```yaml
name: samtools
channels:
  - conda-forge
  - bioconda
dependencies:
  - samtools=1.17
```

**envs/subread.yaml** - For gene counting:
```yaml
name: subread
channels:
  - conda-forge
  - bioconda
dependencies:
  - subread=2.0.3
```

**envs/multiqc.yaml** - For report generation:
```yaml
name: multiqc
channels:
  - conda-forge
  - bioconda
dependencies:
  - multiqc=1.15
```

### Step 3: Configure Your Analysis

Edit the configuration section at the top of the `Snakefile`:

```python
# Reference genome files - Update these paths!
REF_FASTA = "/path/to/your/reference_genome.fa"
REF_GTF = "/path/to/your/gene_annotations.gtf"

# Update your sample information
SAMPLES = {
    "sample1": {"R1": "data/sample1_R1.fastq.gz", "R2": "data/sample1_R2.fastq.gz"},
    "sample2": {"R1": "data/sample2_R1.fastq.gz", "R2": "data/sample2_R2.fastq.gz"},
    # Add all your samples here...
}
```

### Step 4: Prepare Your Data

Place your FASTQ files in the `data/` directory. The pipeline expects paired-end reads with consistent naming:
```
data/
├── sample1_R1.fastq.gz  # Forward reads
├── sample1_R2.fastq.gz  # Reverse reads
├── sample2_R1.fastq.gz
├── sample2_R2.fastq.gz
└── ...
```

## 🏃‍♂️ Running the Pipeline

### Test Run (Dry Run)
Before committing to the full analysis, let's do a practice run to check everything is set up correctly:

```bash
snakemake --dry-run --use-conda
```

This is like doing a dress rehearsal—Snakemake will tell you what it *would* do without actually doing it.

### Full Pipeline Execution

```bash
# Run with 8 CPU cores (adjust based on your computer)
snakemake --use-conda --cores 8
```

**What happens during execution:**
1. Snakemake automatically creates conda environments (first run takes longer)
2. Downloads and installs required software
3. Processes samples in parallel when possible
4. Creates detailed logs for each step
5. Generates comprehensive quality reports

### Monitor Progress

While running, you'll see output like:
```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       build_star_index
        24      fastqc_raw
        24      trim_reads
        ...
```

## 📊 Understanding Your Results

Once complete, explore the `results/` directory:

```
results/
├── qc/
│   ├── raw/              # Quality reports for original data
│   └── trimmed/          # Quality reports after cleaning
├── trimmed/              # Cleaned sequencing reads
├── star_index/           # Genome index files
├── aligned/              # Aligned reads (BAM files)
├── counts/               # Gene expression counts
│   ├── count_matrix.tsv  # 📈 Main results table
│   └── gene_counts.txt   # Detailed count information
└── multiqc/              # 🎨 Beautiful summary report
    └── multiqc_report.html
```

### Key Results Files

**🎯 `count_matrix.tsv`** - Your main results file containing gene expression counts for all samples. Each row is a gene, each column is a sample.

**📊 `multiqc_report.html`** - Open this in your web browser for a comprehensive, interactive quality report with beautiful visualizations.

## 🔍 Quality Control Checkpoints

Good science requires checking your work! Here's what to look for:

### In the MultiQC Report:
- **Sequence quality**: Most bases should have quality scores >30
- **Adapter content**: Should be removed after trimming
- **Alignment rates**: Typically >80% for good samples
- **Gene detection**: Thousands of genes should be detected per sample

### Red Flags to Watch For:
- ⚠️ Very low alignment rates (<50%)
- ⚠️ High duplication levels (>50%)
- ⚠️ Unexpected adapter contamination
- ⚠️ Biased coverage patterns

## 🛠 Troubleshooting Common Issues

### "No such file or directory" errors
- Check that your FASTQ file paths in the `SAMPLES` dictionary are correct
- Ensure reference genome files exist and paths are absolute

### Memory errors during STAR alignment
- Reduce the number of parallel jobs: `snakemake --cores 4`
- Increase memory allocation in the STAR rules

### Conda environment issues
```bash
# Clean conda cache and retry
conda clean --all
snakemake --use-conda --cores 8
```

## 🎨 Customizing Your Analysis

### Adding More Samples
Simply update the `SAMPLES` dictionary in the Snakefile:
```python
SAMPLES = {
    "new_sample": {"R1": "data/new_sample_R1.fastq.gz", "R2": "data/new_sample_R2.fastq.gz"},
    # ... existing samples
}
```

### Changing Parameters
Adjust settings in the configuration section:
```python
# Make trimming more or less stringent
FC_STRANDED = 2  # Change strandedness if needed
READ_LENGTH = 100  # Adjust for your sequencing protocol
```

### Running Specific Steps
```bash
# Only run quality control
snakemake --use-conda results/multiqc/multiqc_report.html

# Only process specific samples
snakemake --use-conda results/aligned/sample1.bam
```

## 📈 Next Steps: Downstream Analysis

With your count matrix ready, you can move on to:
- **Differential expression analysis** (DESeq2, edgeR)
- **Gene pathway analysis** (GSEA, GO enrichment)
- **Visualization** (heatmaps, volcano plots)
- **Machine learning** (clustering, classification)

## 🤝 Getting Help

### Documentation Resources
- [Snakemake Tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [RNA-seq Analysis Best Practices](https://www.nature.com/articles/s41576-019-0150-2)
- [MultiQC Documentation](https://multiqc.info/)

### Community Support
- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Bioconductor Support](https://support.bioconductor.org/)
- [SEQanswers Forums](http://seqanswers.com/)


## 🎉 Congratulations!

You've just set up and run a professional-grade RNA-seq analysis pipeline! This workflow is:
- ✅ **Reproducible**: Anyone can run it with the same results
- ✅ **Scalable**: Works with dozens or hundreds of samples
- ✅ **Documented**: Every step is logged and traceable
- ✅ **Publication-ready**: Generates figures and tables for papers

Remember, bioinformatics is like learning a new language—it takes practice, but each analysis teaches you something new. Don't be discouraged if everything doesn't work perfectly the first time. Even experienced bioinformaticians spend time troubleshooting!

Happy analyzing! 🧬✨



*Found this helpful? Please ⭐ star this repository and share with fellow researchers!*

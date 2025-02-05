# Genome_Alignment_Pipeline

This repository provides a step-by-step guide for genome alignment using breast cancer and non-tumor breast cell line datasets. It includes quality control, read alignment, post-processing, and visualization using various bioinformatics tools.

---

## Table of Contents

1. [Server Information](#1-server-information)  
2. [Dataset](#2-dataset)  
3. [Quality Visualization](#3-quality-visualization)  
4. [Quality Filtering and Trimming](#4-quality-filtering-and-trimming)  
5. [Genome Alignment](#5-genome-alignment)  
6. [Sort and Convert to BAM Format](#6-sort-and-convert-to-bam-format)  
7. [Mark Duplicates](#7-mark-duplicates)  
8. [Recalibrate Base Quality Scores](#8-recalibrate-base-quality-scores)  
9. [Index BAM File](#9-index-bam-file)
10. [Visualization using IGV](#10-visualization-using-igv)  
11. [Exercise](#11-exercise)
    - [11.1 DNAseq Data](#11-1-dnaseq-data)
    - [11.2 QC Raw Reads](#11-2-qc-raw-reads)
    - [11.3 Align Reads using BWA MEM](#11-3-align-reads-using-bwa-mem)
    - [11.4 Sort and Convert to BAM Format](#11-4-sort-and-convert-to-bam-format)
    - [11.5 Mark Duplicates](#11-5-mark-duplicates)
    - [11.6 Recalibrate Base Quality Scores](#11-6-recalibrate-base-quality-scores)
    - [11.7 Index BAM File](11-7-index-bam-file)
13. [Useful Links](#12-useful-links)

---

## 1. Server Information

### Prerequisites
- **Tools:** Install  [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html), [FileZilla](https://filezilla-project.org/download.php), [Bandage](http://rrwick.github.io/Bandage/), [Mauve](https://darlinglab.org/mauve/mauve.html).  
- **Server Details:**  
  - **IP Address:** `168.105.161.70`  
  - **Port:** `22`  
  - **Access:** Requires JABSOM or UH network (use VPN for remote access).
  - Note: With PuTTY and FileZilla you can connect to server.
 
    ![PuTTY](putty_image.PNG "PuTTY") ![FileZilla](FileZilla.png "FileZilla")

### Security Practices
- Avoid multiple failed login attempts to prevent account locking.  
- Use strong passwords or SSH keys.  
- Log out after completing tasks.

### Installing Miniconda
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### Installing software
```bash
# conda install mamba
# mamba install sra-tools fastqc trimmomatic multiqc curl spades quast
mamba install bwa samtools picard gatk4
```

---

## 2. Dataset

```
fastq-dump --split-files -X 100000 SRR097848
```

---


## 3. Quality Visualization

```
mkdir qc
fastqc *.fastq -o qc/
```

---


## 4. Quality Filtering and Trimming

```
curl -OL https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa > adapters.fa
trimmomatic PE SRR097848_1.fastq SRR097848_2.fastq trimmed_1.fastq unpaired_1.fastq trimmed_2.fastq unpaired_2.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:20
```
---



## 5. Genome Alignment

```
bwa index /home/bqhs/hg38/BWAIndex/version0.6.0/genome.fa
bwa mem -R '@RG\tID:SRR097848\tSM:SRR097848\tPL:ILLUMINA\tLB:SRR097848' /home/bqhs/hg38/BWAIndex/version0.6.0/genome.fa trimmed_1.fastq trimmed_2.fastq > SRR097848_raw.sam
```

---



## 6. Sort and Convert to BAM Format

```
samtools sort SRR097848_raw.sam > SRR097848_sort.bam
```


---

## 7. Mark Duplicates  

### Why is it important to mark duplicates in sequencing data?
- **Removes PCR & Optical Duplicates** – Prevents biases from artificially amplified reads.
- **Prevents Overestimation of Read Depth** – Avoids false high coverage that misleads analysis.
- **Improves Variant Calling Accuracy** – Helps distinguish real mutations from PCR artifacts.
- **Reduces Mapping Bias** – Prevents overrepresentation of repetitive/highly expressed regions.

### Why "Mark" Instead of "Remove"? 
- Allows tools (e.g., GATK) to ignore duplicates without losing raw data.
    
```
picard MarkDuplicates -Xmx50g I=SRR097848_sort.bam O=SRR097848_dedup.bam M=SRR097848_dedup.txt
picard CollectAlignmentSummaryMetrics -Xmx50g INPUT=SRR097848_dedup.bam OUTPUT=SRR097848_aln_metrics.txt REFERENCE_SEQUENCE=/home/bqhs/hg38/genome.fa
```

---

## 8. Recalibrate Base Quality Scores  

```
gatk BaseRecalibrator -R /home/bqhs/hg38/genome.fa -I SRR097848_dedup.bam --known-sites /home/bqhs/hg38/dbsnp_146.hg38.vcf.gz --known-sites /home/bqhs/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O recal.table
gatk ApplyBQSR -R /home/bqhs/hg38/genome.fa -I SRR097848_dedup.bam --bqsr-recal-file recal.table -O SRR097848_FINAL.bam
```

---

## 9. Index BAM File

```
samtools index SRR097848_FINAL.bam
```

---

## 10. Visualization using IGV  

Load indexed BAM files in genome browsers IGV: IGV Software

Example Coordinates:

chr1:125,175,342-125,184,719

chr11:47,599,786-47,599,865

---

## 11. Exercise


### 11.1 DNAseq Data

```
mkdir -p breast_cancer
cd breast_cancer
fastq-dump --split-files -X 100000 SRR097849

```
---

### 11.2 QC Raw Reads
```
mkdir -p qc
fastqc *.fastq -o qc/
curl -OL https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa > adapters.fa
trimmomatic PE SRR097849_1.fastq SRR097849_2.fastq trimmed_1.fastq unpaired_1.fastq trimmed_2.fastq unpaired_2.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:20
```
---

### 11.3 Align Reads using BWA MEM
```
bwa mem -t 10 -R '@RG\tID:SRR097849\tSM:SRR097849\tPL:ILLUMINA\tLB:SRbwR097849' /home/bqhs/hg38/BWAIndex/version0.6.0/genome.fa trimmed_1.fastq trimmed_2.fastq > SRR097849_raw.sam
```
---

### 11.4 Sort and Convert to BAM Format

```
samtools sort SRR097849_raw.sam > SRR097849_sort.bam
```
---

### 11.5 Mark Duplicates
```
picard MarkDuplicates -Xmx50g I=SRR097849_sort.bam O=SRR097849_dedup.bam M=SRR097849_dedup.txt
picard CollectAlignmentSummaryMetrics -Xmx50g INPUT=SRR097849_dedup.bam OUTPUT=SRR097849_aln_metrics.txt REFERENCE_SEQUENCE=/home/bqhs/hg38/genome.fa
```
---

### 11.6 Recalibrate Base Quality Scores
```
gatk BaseRecalibrator -R /home/bqhs/hg38/genome.fa -I SRR097849_dedup.bam --known-sites /home/bqhs/hg38/dbsnp_146.hg38.vcf.gz --known-sites /home/bqhs/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O recal.table
gatk ApplyBQSR -R /home/bqhs/hg38/genome.fa -I SRR097849_dedup.bam --bqsr-recal-file recal.table -O SRR097849_FINAL.bam
```
---

### 11.7 Index BAM File
```
samtools index SRR097849_FINAL.bam
```
---

---

## 12. Useful Links


---



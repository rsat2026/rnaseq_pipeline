# RNA-seq Workflow (STAR + featureCounts)

This package provides a beginner-friendly RNA-seq pipeline.

It supports:
- **Paired-end** FASTQ.gz (`*_R1*` + `*_R2*` or `_1/_2`)
- **Single-end** FASTQ.gz (no mate)

Outputs:
- per-sample QC (FastQC)
- trimming reports (fastp)
- aligned BAMs (STAR + samtools index)
- gene counts (featureCounts)
- MultiQC report
- `SUMMARY.tsv` per sample

## Cross-platform notes

- **Linux/Ubuntu:** Works directly using `bash`.
- **macOS:** Works, but you must run with **Bash 4+** (macOS default bash is often 3.2).
  - Install brew bash: `brew install bash gawk`
  - Run script using brew bash path (example): `/opt/homebrew/bin/bash rnaseq_pipeline.sh ...`
- **Windows:** Use **WSL2 + Ubuntu**, and run everything inside the Ubuntu terminal.

---

## Minimum PC requirements (personal computer)

These depend heavily on genome size and sample count.

- **CPU:** 4 cores minimum (8+ recommended)
- **RAM:**
  - small genomes (yeast/bacteria): **8–16 GB** often OK
  - human/mouse: **32 GB recommended** for smoother STAR indexing/mapping
- **Disk:** **50–150+ GB free** (FASTQs + BAMs + STAR outputs + QC)
- **Internet:** needed to download tools and reference files

---

## Project structure

Recommended layout (STAR index directly under the project folder):
# STAR Genome Index Build Instructions (RNA-seq Workflow)

This document explains how to build the STAR genome index required by the RNA-seq pipeline.

The STAR index must be created before running:  rnaseq_pipeline.sh

---
## Overview of STAR index builder##

Two STAR index builder scripts are provided:

| Script | When to use |
|-------|------------|
| `star_index_build.sh` | Recommended default (faster, better performance if you have enough RAM) |
| `star_index_build_lowram.sh` | Use this if STAR fails due to memory limits |

Both scripts:

- Download tested GENCODE v44 human references
- Create the required project structure
- Build the index into: reference/STAR_index

This matches the pipeline defaults.

---
## Step 0 — Activate environment

Before building the index, activate the RNA-seq environment: mamba activate rnaseq

If you do not have it yet:  

mamba env create -f env.yml
mamba activate rnaseq
---

## Step 1 — Choose your STAR build method

### Recommended (standard)

Use this if your computer has:

- 32 GB RAM (ideal)
- OR 16 GB RAM (often works)

Run:  bash star_index_build.sh -p my_rnaseq_project -t 8 -o 100

---

### Low-RAM version

Use this if you see errors like: EXITING because of FATAL ERROR: not enough memory


Run:  bash star_index_build_lowram.sh -p my_rnaseq_project -t 4 -o 100


This version:

- reduces STAR index memory usage
- may take longer to build
- may slightly reduce alignment performance
- but works better on personal laptops

---

## Step 2 — Understand sjdbOverhang

`-o` controls:   sjdbOverhang = read_length - 1


Examples:

| Read length | Use |
|------------|-----|
| 101 bp | `-o 100` |
| 151 bp | `-o 150` |
| unknown | use `100` (safe default) |

---

## Step 3 — What files will be downloaded

The scripts automatically download tested references:

gencode.v44.primary_assembly.annotation.gtf
GRCh38.primary_assembly.genome.fa


These are placed in:  my_rnaseq_project/reference/

---

## Step 4 — Expected project structure after build


my_rnaseq_project/
fastq/
reference/
GRCh38.primary_assembly.genome.fa
gencode.v44.primary_assembly.annotation.gtf
STAR_index/
output/

---

## Step 5 — Run the RNA-seq pipeline

After STAR index finishes:


cd my_rnaseq_project

bash rnaseq_pipeline.sh
-i fastq
-o output
-g reference/STAR_index
-a reference/gencode.v44.primary_assembly.annotation.gtf
-t 4
--stranded 0


---

## Troubleshooting

### STAR build fails due to RAM

Try: bash star_index_build_lowram.sh -t 2

Close other programs.

STAR indexing for human genomes can require large memory.

---

### STAR command not found

You forgot to activate environment:

mamba activate rnaseq

---

### Download fails

Check internet connection.

You can manually download the two files and place them in:

reference/

Then rerun the script.

---

## Notes for macOS users

macOS default bash is often outdated.

Install modern bash:

brew install bash gawk

Run scripts using:


/opt/homebrew/bin/bash star_index_build.sh ...
---

## Final recommendation

- Use the standard build if possible
- Use low-RAM build only if necessary

Both produce indexes compatible with the RNA-seq pipeline.
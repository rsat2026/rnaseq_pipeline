#!/usr/bin/env bash
# star_index_build.sh
# Download GENCODE v44 primary assembly FASTA/GTF (GRCh38) and build STAR index.
# Tested file names:
#   gencode.v44.primary_assembly.annotation.gtf
#   GRCh38.primary_assembly.genome.fa
#
# Usage:
#   bash star_index_build.sh -p my_rnaseq_project -t 8 -m 16 -o 100
#
# Options:
#   -p  Project directory (default: my_rnaseq_project)
#   -t  Threads for STAR (default: 8)
#   -m  RAM limit in GB for STAR genomeGenerate (default: 16)
#   -o  sjdbOverhang (read_length-1). For 101bp reads use 100 (default: 100)
#
# Requirements:
#   STAR, wget, gunzip (gzip), and a working conda/mamba environment is recommended.

set -euo pipefail

PROJECT_DIR="my_rnaseq_project"
THREADS=8
RAM_GB=16
SJDB_OVERHANG=100

usage() {
  cat <<'USAGE'
Usage:
  bash star_index_build.sh [-p PROJECT_DIR] [-t THREADS] [-m RAM_GB] [-o SJDB_OVERHANG]

Examples:
  bash star_index_build.sh -p my_rnaseq_project -t 8 -m 16 -o 100
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p) PROJECT_DIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -m) RAM_GB="$2"; shift 2 ;;
    -o) SJDB_OVERHANG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if ! command -v STAR >/dev/null 2>&1; then
  echo "ERROR: STAR not found in PATH. Activate your conda env first (mamba activate rnaseq)." >&2
  exit 1
fi

mkdir -p "$PROJECT_DIR"/{fastq,reference,output}
cd "$PROJECT_DIR/reference"

GTF_GZ="gencode.v44.primary_assembly.annotation.gtf.gz"
FA_GZ="GRCh38.primary_assembly.genome.fa.gz"
GTF="gencode.v44.primary_assembly.annotation.gtf"
FA="GRCh38.primary_assembly.genome.fa"
STAR_DIR="STAR_index"

GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
FA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"

echo "[INFO] Project: $PROJECT_DIR"
echo "[INFO] Threads: $THREADS"
echo "[INFO] RAM limit: ${RAM_GB} GB"
echo "[INFO] sjdbOverhang: $SJDB_OVERHANG"
echo "[INFO] Downloading reference files into: $PROJECT_DIR/reference"

if [[ ! -f "$GTF_GZ" ]]; then
  echo "[INFO] Downloading GTF..."
  wget -O "$GTF_GZ" "$GTF_URL"
else
  echo "[INFO] GTF.gz exists, skipping download."
fi

if [[ ! -f "$FA_GZ" ]]; then
  echo "[INFO] Downloading FASTA..."
  wget -O "$FA_GZ" "$FA_URL"
else
  echo "[INFO] FASTA.gz exists, skipping download."
fi

if [[ ! -f "$GTF" ]]; then
  echo "[INFO] Unzipping GTF..."
  gunzip -c "$GTF_GZ" > "$GTF"
else
  echo "[INFO] GTF exists, skipping unzip."
fi

if [[ ! -f "$FA" ]]; then
  echo "[INFO] Unzipping FASTA..."
  gunzip -c "$FA_GZ" > "$FA"
else
  echo "[INFO] FASTA exists, skipping unzip."
fi

mkdir -p "$STAR_DIR"

# Convert GB to bytes for STAR parameter
LIMIT_RAM_BYTES=$((RAM_GB * 1024 * 1024 * 1024))

echo "[INFO] Building STAR index..."
STAR \
  --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_DIR" \
  --genomeFastaFiles "$FA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --genomeSAindexNbases 14 \
  --limitGenomeGenerateRAM "$LIMIT_RAM_BYTES"

echo "[DONE] STAR index built at: $PROJECT_DIR/reference/$STAR_DIR"
echo "[NEXT] Run the pipeline like:"
echo "  bash rnaseq_pipeline.sh -i fastq -o output -g reference/STAR_index -a reference/$GTF -t $THREADS --stranded 0"
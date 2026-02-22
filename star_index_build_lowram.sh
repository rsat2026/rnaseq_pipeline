#!/usr/bin/env bash
# star_index_build_lowram.sh
# Low-RAM STAR genome index builder (GRCh38 + GENCODE v44 primary assembly)
#
# Builds STAR index into:
#   <project>/reference/STAR_index
#
# Uses low-RAM tuning:
#   --genomeSAindexNbases 12
#   --genomeSAsparseD 2
#   --genomeChrBinNbits 18
#
# Also includes RNA-seq recommended annotation:
#   --sjdbGTFfile + --sjdbOverhang
#
# Usage:
#   bash star_index_build_lowram.sh -p my_rnaseq_project -t 4 -o 100
#
# Notes:
# - sjdbOverhang should be (read_length - 1). For 101bp reads use 100; for 151bp reads use 150.
# - Low-RAM settings may reduce performance slightly, but can help on personal PCs with limited memory.

set -euo pipefail

PROJECT_DIR="my_rnaseq_project"
THREADS=4
SJDB_OVERHANG=100

usage() {
  cat <<'USAGE'
Usage:
  bash star_index_build_lowram.sh [-p PROJECT_DIR] [-t THREADS] [-o SJDB_OVERHANG]

Options:
  -p  Project directory (default: my_rnaseq_project)
  -t  Threads (default: 4)  # lower threads often helps reduce overall memory pressure
  -o  sjdbOverhang = read_length - 1 (default: 100)

Examples:
  bash star_index_build_lowram.sh -p my_rnaseq_project -t 4 -o 100
  bash star_index_build_lowram.sh -p my_rnaseq_project -t 4 -o 150
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p) PROJECT_DIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -o) SJDB_OVERHANG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

need_cmd() {
  local c="$1"
  if ! command -v "$c" >/dev/null 2>&1; then
    echo "ERROR: Required command not found: $c" >&2
    exit 1
  fi
}

need_cmd STAR
need_cmd wget
need_cmd gunzip
need_cmd tee

mkdir -p "$PROJECT_DIR"/{fastq,reference,output}
cd "$PROJECT_DIR/reference"

GTF_GZ="gencode.v44.primary_assembly.annotation.gtf.gz"
FA_GZ="GRCh38.primary_assembly.genome.fa.gz"
GTF="gencode.v44.primary_assembly.annotation.gtf"
FA="GRCh38.primary_assembly.genome.fa"
STAR_DIR="STAR_index"

GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
FA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"

LOGFILE="../output/star_genomeGenerate.lowram.log"

echo "[INFO] Project directory: $PROJECT_DIR"
echo "[INFO] Reference dir: $(pwd)"
echo "[INFO] Threads: $THREADS"
echo "[INFO] sjdbOverhang: $SJDB_OVERHANG"
echo "[INFO] STAR index dir: $PROJECT_DIR/reference/$STAR_DIR"
echo "[INFO] Log: $LOGFILE"
echo

# Download tested references (skip if already present)
if [[ ! -f "$GTF_GZ" ]]; then
  echo "[INFO] Downloading GTF..."
  wget -O "$GTF_GZ" "$GTF_URL"
else
  echo "[INFO] Found $GTF_GZ (skip download)"
fi

if [[ ! -f "$FA_GZ" ]]; then
  echo "[INFO] Downloading FASTA..."
  wget -O "$FA_GZ" "$FA_URL"
else
  echo "[INFO] Found $FA_GZ (skip download)"
fi

# Unzip (skip if already unzipped)
if [[ ! -f "$GTF" ]]; then
  echo "[INFO] Unzipping GTF..."
  gunzip -c "$GTF_GZ" > "$GTF"
else
  echo "[INFO] Found $GTF (skip unzip)"
fi

if [[ ! -f "$FA" ]]; then
  echo "[INFO] Unzipping FASTA..."
  gunzip -c "$FA_GZ" > "$FA"
else
  echo "[INFO] Found $FA (skip unzip)"
fi

mkdir -p "$STAR_DIR"
mkdir -p ../output

echo "[INFO] Building STAR index with LOW-RAM settings..."
echo "[INFO] If this still fails due to memory, reduce threads (-t 2) and close other programs."
echo

# Low-RAM STAR genomeGenerate (your tuning) + annotation for RNA-seq
STAR \
  --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_DIR" \
  --genomeFastaFiles "$FA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --genomeSAindexNbases 12 \
  --genomeSAsparseD 2 \
  --genomeChrBinNbits 18 \
  2>&1 | tee "$LOGFILE"

echo
echo "[DONE] STAR index built at: $PROJECT_DIR/reference/$STAR_DIR"
echo "[DONE] Log saved to: $PROJECT_DIR/output/star_genomeGenerate.lowram.log"
echo
echo "[NEXT] Run the pipeline:"
echo "  cd $PROJECT_DIR"
echo "  bash rnaseq_pipeline.sh -i fastq -o output -g reference/STAR_index -a reference/$GTF -t $THREADS --stranded 0"
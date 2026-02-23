#!/usr/bin/env bash
# rnaseq_pipeline.sh
# Production-ready RNA-seq pipeline for SE/PE FASTQ.gz inputs
#
# Usage:
#   bash rnaseq_pipeline.sh -i <INPUT_DIR> -o <OUTPUT_DIR> -g <GENOME_DIR> -a <GTF> [-t THREADS] [--single-end] [--paired-end] [--stranded 0/1/2] [--force] [--dry-run]
#
# Notes:
# - If neither --single-end nor --paired-end is provided, the script auto-detects per sample.
# - Paired-end patterns supported: *_R1*.fastq.gz + *_R2*.fastq.gz or *_1*.fastq.gz + *_2*.fastq.gz
# - Single-end: *.fastq.gz with no R2 partner
# - Requires: fastqc, multiqc, fastp, STAR, samtools, featureCounts

set -euo pipefail

SCRIPT_NAME="rnaseq_pipeline.sh"

# -----------------------------
# Defaults
# -----------------------------
THREADS=4
READ_MODE="auto"   # auto|se|pe
STRANDED=0
FORCE=0
DRY_RUN=0

INPUT_DIR=""
OUTPUT_DIR=""
GENOME_DIR=""
GTF=""

# -----------------------------
# Helpers
# -----------------------------
usage() {
  cat <<'USAGE'
Usage:
  bash rnaseq_pipeline.sh -i <INPUT_DIR> -o <OUTPUT_DIR> -g <GENOME_DIR> -a <GTF> [-t THREADS] [--single-end] [--paired-end] [--stranded 0/1/2] [--force] [--dry-run]

Options:
  -i, --input       Input directory containing FASTQ.gz files
  -o, --output      Output directory
  -g, --genome      STAR genome index directory
  -a, --gtf         GTF annotation file
  -t, --threads     Threads (default: 4)
  --single-end      Force single-end mode
  --paired-end      Force paired-end mode
  --stranded        Strandedness for featureCounts: 0 (unstranded), 1 (stranded), 2 (reversely stranded). Default: 0
  --force           Overwrite/redo existing outputs
  --dry-run         Print actions without executing
  -h, --help        Show this help
USAGE
}

log() {
  local msg="$1"
  local ts
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[$ts] $msg"
  else
    echo "[$ts] $msg" | tee -a "$OUTPUT_DIR/logs/pipeline.log" >/dev/null
  fi
}

run_cmd() {
  local cmd="$1"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] $cmd"
    return 0
  fi
  eval "$cmd"
}

quote_args() {
  local out=()
  local a
  for a in "$@"; do
    out+=("$(printf '%q' "$a")")
  done
  echo "${out[*]}"
}

on_error() {
  local lineno="$1"
  local cmd="$2"
  echo "[ERROR] Line $lineno: $cmd" >&2
}
trap 'on_error ${LINENO} "$BASH_COMMAND"' ERR

# -----------------------------
# Arg parsing
# -----------------------------
if [[ "$#" -eq 0 ]]; then
  usage
  exit 1
fi

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT_DIR="$2"; shift 2 ;;
    -o|--output)
      OUTPUT_DIR="$2"; shift 2 ;;
    -g|--genome)
      GENOME_DIR="$2"; shift 2 ;;
    -a|--gtf)
      GTF="$2"; shift 2 ;;
    -t|--threads)
      THREADS="$2"; shift 2 ;;
    --single-end)
      READ_MODE="se"; shift ;;
    --paired-end)
      READ_MODE="pe"; shift ;;
    --stranded)
      STRANDED="$2"; shift 2 ;;
    --force)
      FORCE=1; shift ;;
    --dry-run)
      DRY_RUN=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$GENOME_DIR" || -z "$GTF" ]]; then
  echo "Error: -i, -o, -g, and -a are required." >&2
  usage
  exit 1
fi

if [[ "$STRANDED" != "0" && "$STRANDED" != "1" && "$STRANDED" != "2" ]]; then
  echo "Error: --stranded must be 0, 1, or 2." >&2
  exit 1
fi

# -----------------------------
# Preflight checks
# -----------------------------
if [[ "$DRY_RUN" -eq 0 ]]; then
  for tool in fastqc multiqc fastp STAR samtools featureCounts; do
    if ! command -v "$tool" >/dev/null 2>&1; then
      echo "Error: Required tool not found in PATH: $tool" >&2
      exit 1
    fi
  done
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Input dir not found: $INPUT_DIR" >&2
  exit 1
fi

if [[ ! -d "$GENOME_DIR" ]]; then
  echo "Error: STAR genome dir not found: $GENOME_DIR" >&2
  exit 1
fi

if [[ ! -f "$GTF" ]]; then
  echo "Error: GTF not found: $GTF" >&2
  exit 1
fi

# Check STAR index files
STAR_REQ=("Genome" "SA" "SAindex" "chrLength.txt" "chrName.txt" "chrStart.txt" "genomeParameters.txt")
for f in "${STAR_REQ[@]}"; do
  if [[ ! -f "$GENOME_DIR/$f" ]]; then
    echo "Error: STAR index missing: $GENOME_DIR/$f" >&2
    exit 1
  fi
done

# Verify input FASTQ.gz exists
shopt -s nullglob
FASTQS=("$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz)
shopt -u nullglob
if [[ "${#FASTQS[@]}" -eq 0 ]]; then
  echo "Error: No FASTQ.gz files found in $INPUT_DIR" >&2
  exit 1
fi

# -----------------------------
# Output structure
# -----------------------------
run_cmd "mkdir -p \"$OUTPUT_DIR\""
run_cmd "mkdir -p \"$OUTPUT_DIR/logs\" \"$OUTPUT_DIR/qc_raw\" \"$OUTPUT_DIR/qc_trimmed\" \"$OUTPUT_DIR/trimmed\" \"$OUTPUT_DIR/star\" \"$OUTPUT_DIR/bam\" \"$OUTPUT_DIR/counts\" \"$OUTPUT_DIR/multiqc\""

log "Starting RNA-seq pipeline"
log "Input: $INPUT_DIR"
log "Output: $OUTPUT_DIR"
log "Genome: $GENOME_DIR"
log "GTF: $GTF"
log "Threads: $THREADS"
log "Read mode: $READ_MODE"
log "Strandedness: $STRANDED"
log "Force: $FORCE"
log "Dry-run: $DRY_RUN"

# -----------------------------
# Sample discovery
# -----------------------------
if [[ ${BASH_VERSINFO[0]} -lt 4 ]]; then
  echo "Error: Bash 4+ required for associative arrays." >&2
  exit 1
fi

declare -A R1_MAP
declare -A R2_MAP
declare -A SE_MAP

detect_fastqs() {
  local f base sample
  local re_r1='^(.+)(_R1|_1)([^/]*)\.(f(ast)?q)\.gz$'
  local re_r2='^(.+)(_R2|_2)([^/]*)\.(f(ast)?q)\.gz$'

  for f in "${FASTQS[@]}"; do
    base=$(basename "$f")
    if [[ "$base" =~ $re_r1 ]]; then
      sample="${BASH_REMATCH[1]}"
      R1_MAP["$sample"]="$f"
    elif [[ "$base" =~ $re_r2 ]]; then
      sample="${BASH_REMATCH[1]}"
      R2_MAP["$sample"]="$f"
    else
      # Single-end candidate
      sample="$base"
      sample="${sample%.fastq.gz}"
      sample="${sample%.fq.gz}"
      SE_MAP["$sample"]="$f"
    fi
  done
}

detect_fastqs

# Build sample list
SAMPLES=()

declare -A SAMPLE_SET
for s in "${!R1_MAP[@]}" "${!R2_MAP[@]}" "${!SE_MAP[@]}"; do
  if [[ -n "$s" ]]; then
    SAMPLE_SET["$s"]=1
  fi
done

for s in "${!SAMPLE_SET[@]}"; do
  SAMPLES+=("$s")
done

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
  echo "Error: No samples detected." >&2
  exit 1
fi

# Sort samples for deterministic order
IFS=$'\n' SAMPLES=($(sort <<<"${SAMPLES[*]}"))
unset IFS

# Summary file
SUMMARY_TSV="$OUTPUT_DIR/SUMMARY.tsv"
if [[ "$FORCE" -eq 1 || ! -f "$SUMMARY_TSV" ]]; then
  run_cmd "printf 'sample_id\tread_type\traw_reads\ttrimmed_reads\tstar_unique_pct\tstar_multi_pct\tassigned_reads\terror_flags\n' > \"$SUMMARY_TSV\""
fi

# Helper to parse fastp json
parse_fastp_reads() {
  local json="$1"
  if [[ ! -f "$json" ]]; then
    echo -e "\t"
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    python3 - <<'PY' "$json"
import json,sys
p=sys.argv[1]
try:
    data=json.load(open(p))
    before=data.get('summary',{}).get('before_filtering',{}).get('total_reads','')
    after=data.get('summary',{}).get('after_filtering',{}).get('total_reads','')
    print(f"{before}\t{after}")
except Exception:
    print("\t")
PY
  else
    echo -e "\t"
  fi
}

# Helper to parse STAR Log.final.out
parse_star_log() {
  local logf="$1"
  local uniq="" multi=""
  if [[ -f "$logf" ]]; then
    uniq=$(awk -F '|' '/Uniquely mapped reads %/ {gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2}' "$logf" | head -n1)
    multi=$(awk -F '|' '/% of reads mapped to multiple loci/ {gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2}' "$logf" | head -n1)
  fi
  echo -e "${uniq}\t${multi}"
}

# Helper to parse featureCounts assigned reads
parse_fc_assigned() {
  local summary="$1"
  local assigned=""
  if [[ -f "$summary" ]]; then
    assigned=$(awk -F '\t' '$1=="Assigned" {print $2}' "$summary" | head -n1)
  fi
  echo "$assigned"
}

# Write summary row
write_summary() {
  local sample="$1"
  local read_type="$2"
  local error_flags="$3"
  local fastp_json="$4"
  local star_log="$5"
  local counts_summary="$6"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] SUMMARY: $sample\t$read_type\t...\t$error_flags"
    return 0
  fi

  local raw_trimmed raw_reads trimmed_reads
  raw_trimmed=$(parse_fastp_reads "$fastp_json")
  raw_reads=$(echo -e "$raw_trimmed" | awk -F '\t' '{print $1}')
  trimmed_reads=$(echo -e "$raw_trimmed" | awk -F '\t' '{print $2}')

  local star_stats star_unique star_multi
  star_stats=$(parse_star_log "$star_log")
  star_unique=$(echo -e "$star_stats" | awk -F '\t' '{print $1}')
  star_multi=$(echo -e "$star_stats" | awk -F '\t' '{print $2}')

  local assigned
  assigned=$(parse_fc_assigned "$counts_summary")

  echo -e "$sample\t$read_type\t$raw_reads\t$trimmed_reads\t$star_unique\t$star_multi\t$assigned\t$error_flags" >> "$SUMMARY_TSV"
}

# -----------------------------
# Per-sample processing
# -----------------------------
BAM_LIST=()

process_sample() {
  local sample="$1"
  local error_flags=""
  local read_type=""
  local r1="${R1_MAP[$sample]:-}"
  local r2="${R2_MAP[$sample]:-}"
  local se="${SE_MAP[$sample]:-}"

  # Determine read type
  if [[ "$READ_MODE" == "se" ]]; then
    if [[ -n "$se" ]]; then
      read_type="SE"
    elif [[ -n "$r1" && -z "$r2" ]]; then
      read_type="SE"
      se="$r1"
      r1=""
    else
      error_flags="FORCED_SE_MISSING"
      log "[$sample] Error: --single-end but no SE FASTQ found"
      write_summary "$sample" "SE" "$error_flags" "" "" ""
      return 1
    fi
  elif [[ "$READ_MODE" == "pe" ]]; then
    if [[ -n "$r1" && -n "$r2" ]]; then
      read_type="PE"
    else
      error_flags="FORCED_PE_MISSING_MATE"
      log "[$sample] Error: --paired-end but missing mate"
      write_summary "$sample" "PE" "$error_flags" "" "" ""
      return 1
    fi
  else
    # Auto-detect per sample
    if [[ -n "$r1" || -n "$r2" ]]; then
      if [[ -n "$r1" && -n "$r2" ]]; then
        read_type="PE"
      else
        error_flags="MISSING_MATE"
        log "[$sample] Missing mate for paired-end pattern"
        write_summary "$sample" "PE" "$error_flags" "" "" ""
        return 1
      fi
    elif [[ -n "$se" ]]; then
      read_type="SE"
    else
      error_flags="NO_INPUT"
      log "[$sample] No input FASTQ found"
      write_summary "$sample" "" "$error_flags" "" "" ""
      return 1
    fi
  fi

  log "[$sample] Processing ($read_type)"

  # Paths
  local qc_raw="$OUTPUT_DIR/qc_raw"
  local qc_trim="$OUTPUT_DIR/qc_trimmed"
  local trimmed="$OUTPUT_DIR/trimmed"
  local star_dir="$OUTPUT_DIR/star/$sample"
  local bam_dir="$OUTPUT_DIR/bam"
  local counts_dir="$OUTPUT_DIR/counts"

  local fastp_html="$qc_trim/${sample}.fastp.html"
  local fastp_json="$qc_trim/${sample}.fastp.json"

  local trim_r1="$trimmed/${sample}_R1.trimmed.fastq.gz"
  local trim_r2="$trimmed/${sample}_R2.trimmed.fastq.gz"
  local trim_se="$trimmed/${sample}.trimmed.fastq.gz"

  local star_bam="$star_dir/Aligned.sortedByCoord.out.bam"
  local final_bam="$bam_dir/${sample}.bam"
  local final_bai="$final_bam.bai"
  local star_log="$star_dir/Log.final.out"

  local counts_file="$counts_dir/${sample}.counts.txt"
  local counts_summary="$counts_dir/${sample}.counts.txt.summary"

  # 1) FastQC on raw reads
  if [[ "$read_type" == "PE" ]]; then
    local r1_base r2_base raw_qc1 raw_qc2
    r1_base=$(basename "$r1"); r1_base=${r1_base%.fastq.gz}; r1_base=${r1_base%.fq.gz}
    r2_base=$(basename "$r2"); r2_base=${r2_base%.fastq.gz}; r2_base=${r2_base%.fq.gz}
    raw_qc1="$qc_raw/${r1_base}_fastqc.html"
    raw_qc2="$qc_raw/${r2_base}_fastqc.html"
    if [[ "$FORCE" -eq 1 || ! -f "$raw_qc1" || ! -f "$raw_qc2" ]]; then
      if ! run_cmd "fastqc -t $THREADS -o \"$qc_raw\" \"$r1\" \"$r2\""; then
        error_flags+="FASTQC_RAW_FAIL;"
      fi
    else
      log "[$sample] Skip FastQC raw (exists)"
    fi
  else
    local se_base raw_qc
    se_base=$(basename "$se"); se_base=${se_base%.fastq.gz}; se_base=${se_base%.fq.gz}
    raw_qc="$qc_raw/${se_base}_fastqc.html"
    if [[ "$FORCE" -eq 1 || ! -f "$raw_qc" ]]; then
      if ! run_cmd "fastqc -t $THREADS -o \"$qc_raw\" \"$se\""; then
        error_flags+="FASTQC_RAW_FAIL;"
      fi
    else
      log "[$sample] Skip FastQC raw (exists)"
    fi
  fi

  # 2) fastp trimming
  if [[ "$read_type" == "PE" ]]; then
    if [[ "$FORCE" -eq 1 || ! -f "$trim_r1" || ! -f "$trim_r2" ]]; then
      if ! run_cmd "fastp -i \"$r1\" -I \"$r2\" -o \"$trim_r1\" -O \"$trim_r2\" --detect_adapter_for_pe -q 20 -l 30 -w $THREADS -h \"$fastp_html\" -j \"$fastp_json\""; then
        error_flags+="FASTP_FAIL;"
        write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
        return 1
      fi
    else
      log "[$sample] Skip fastp (exists)"
    fi
    if [[ ! -f "$trim_r1" || ! -f "$trim_r2" ]]; then
      error_flags+="FASTP_OUTPUT_MISSING;"
      write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
      return 1
    fi
  else
    if [[ "$FORCE" -eq 1 || ! -f "$trim_se" ]]; then
      if ! run_cmd "fastp -i \"$se\" -o \"$trim_se\" -q 20 -l 30 -w $THREADS -h \"$fastp_html\" -j \"$fastp_json\""; then
        error_flags+="FASTP_FAIL;"
        write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
        return 1
      fi
    else
      log "[$sample] Skip fastp (exists)"
    fi
    if [[ ! -f "$trim_se" ]]; then
      error_flags+="FASTP_OUTPUT_MISSING;"
      write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
      return 1
    fi
  fi

  # 3) FastQC on trimmed reads
  if [[ "$read_type" == "PE" ]]; then
    local t1_base t2_base trim_qc1 trim_qc2
    t1_base=$(basename "$trim_r1"); t1_base=${t1_base%.fastq.gz}; t1_base=${t1_base%.fq.gz}
    t2_base=$(basename "$trim_r2"); t2_base=${t2_base%.fastq.gz}; t2_base=${t2_base%.fq.gz}
    trim_qc1="$qc_trim/${t1_base}_fastqc.html"
    trim_qc2="$qc_trim/${t2_base}_fastqc.html"
    if [[ "$FORCE" -eq 1 || ! -f "$trim_qc1" || ! -f "$trim_qc2" ]]; then
      if ! run_cmd "fastqc -t $THREADS -o \"$qc_trim\" \"$trim_r1\" \"$trim_r2\""; then
        error_flags+="FASTQC_TRIM_FAIL;"
      fi
    else
      log "[$sample] Skip FastQC trimmed (exists)"
    fi
  else
    local t_base trim_qc
    t_base=$(basename "$trim_se"); t_base=${t_base%.fastq.gz}; t_base=${t_base%.fq.gz}
    trim_qc="$qc_trim/${t_base}_fastqc.html"
    if [[ "$FORCE" -eq 1 || ! -f "$trim_qc" ]]; then
      if ! run_cmd "fastqc -t $THREADS -o \"$qc_trim\" \"$trim_se\""; then
        error_flags+="FASTQC_TRIM_FAIL;"
      fi
    else
      log "[$sample] Skip FastQC trimmed (exists)"
    fi
  fi

  # 4) Align with STAR
  if [[ "$FORCE" -eq 1 || ! -f "$star_bam" ]]; then
    run_cmd "mkdir -p \"$star_dir\""
    if [[ "$read_type" == "PE" ]]; then
      if ! run_cmd "STAR --runThreadN $THREADS --genomeDir \"$GENOME_DIR\" --readFilesIn \"$trim_r1\" \"$trim_r2\" --readFilesCommand zcat --outFileNamePrefix \"$star_dir/\" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000  --outFilterMultimapNmax 20 --outFilterMismatchNoverReadLmax 0.04 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMstrandField intronMotif"; then
        error_flags+="STAR_FAIL;"
        write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
        return 1
      fi
    else
      if ! run_cmd "STAR --runThreadN $THREADS --genomeDir \"$GENOME_DIR\" --readFilesIn \"$trim_se\" --readFilesCommand zcat --outFileNamePrefix \"$star_dir/\" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000  --outFilterMultimapNmax 20 --outFilterMismatchNoverReadLmax 0.04 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMstrandField intronMotif"; then
        error_flags+="STAR_FAIL;"
        write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
        return 1
      fi
    fi
  else
    log "[$sample] Skip STAR (exists)"
  fi

  if [[ ! -f "$star_bam" ]]; then
    error_flags+="STAR_OUTPUT_MISSING;"
    write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
    return 1
  fi

  # 5) samtools index BAM (copy bam to bam/)
  if [[ "$FORCE" -eq 1 || ! -f "$final_bam" ]]; then
    if ! run_cmd "cp \"$star_bam\" \"$final_bam\""; then
      error_flags+="BAM_COPY_FAIL;"
      write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
      return 1
    fi
  else
    log "[$sample] Skip BAM copy (exists)"
  fi

  if [[ ! -f "$final_bam" ]]; then
    error_flags+="BAM_MISSING;"
    write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"
    return 1
  fi

  if [[ "$FORCE" -eq 1 || ! -f "$final_bai" ]]; then
    if ! run_cmd "samtools index -@ $THREADS \"$final_bam\""; then
      error_flags+="BAM_INDEX_FAIL;"
    fi
  else
    log "[$sample] Skip BAM index (exists)"
  fi

  # 6) featureCounts per sample
  if [[ "$FORCE" -eq 1 || ! -f "$counts_file" ]]; then
    local fc_flags=""
    if [[ "$read_type" == "PE" ]]; then
      fc_flags="-p -B -C"
    fi
    if ! run_cmd "featureCounts -T $THREADS -s $STRANDED $fc_flags -a \"$GTF\" -o \"$counts_file\" \"$final_bam\""; then
      error_flags+="FEATURECOUNTS_FAIL;"
    fi
  else
    log "[$sample] Skip featureCounts (exists)"
  fi

  if [[ ! -f "$counts_file" ]]; then
    error_flags+="FEATURECOUNTS_OUTPUT_MISSING;"
  fi

  # Collect for combined counts
  if [[ -f "$final_bam" ]]; then
    BAM_LIST+=("$final_bam")
  fi

  write_summary "$sample" "$read_type" "$error_flags" "$fastp_json" "$star_log" "$counts_summary"

  if [[ -n "$error_flags" ]]; then
    return 1
  fi
  return 0
}

# -----------------------------
# Run pipeline
# -----------------------------
for sample in "${SAMPLES[@]}"; do
  set +e
  process_sample "$sample"
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    log "[$sample] Completed with errors (see SUMMARY.tsv)"
  else
    log "[$sample] Completed"
  fi
done

# -----------------------------
# Combined counts matrix
# -----------------------------
COMBINED_COUNTS="$OUTPUT_DIR/counts/combined_counts.txt"

if [[ "${#BAM_LIST[@]}" -gt 0 ]]; then
  if [[ "$FORCE" -eq 1 || ! -f "$COMBINED_COUNTS" ]]; then
    bam_args=$(quote_args "${BAM_LIST[@]}")

    # Auto-detect PE vs SE from the first BAM
    # If there are any paired reads (-f 1), treat as paired-end and use -p -B -C.
    FC_COMBINED_FLAGS=""
    if command -v samtools >/dev/null 2>&1; then
      pe_reads=$(samtools view -c -f 1 "${BAM_LIST[0]}" 2>/dev/null || echo 0)
      if [[ "${pe_reads:-0}" -gt 0 ]]; then
        FC_COMBINED_FLAGS="-p -B -C"
      fi
    fi

    run_cmd "featureCounts -T $THREADS -s $STRANDED $FC_COMBINED_FLAGS -a \"$GTF\" -o \"$COMBINED_COUNTS\" $bam_args"
  else
    log "Skip combined featureCounts (exists)"
  fi
else
  log "No BAMs available for combined counts"
fi

# -----------------------------
# MultiQC
# -----------------------------
if [[ "$FORCE" -eq 1 || ! -f "$OUTPUT_DIR/multiqc/multiqc_report.html" ]]; then
  run_cmd "multiqc -o \"$OUTPUT_DIR/multiqc\" \"$OUTPUT_DIR/qc_raw\" \"$OUTPUT_DIR/qc_trimmed\" \"$OUTPUT_DIR/star\" \"$OUTPUT_DIR/counts\""
else
  log "Skip MultiQC (exists)"
fi

log "Pipeline complete"

#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/Breast-Cancer-RNAseq-DE-Analysis"

# Test list first. Later we will switch this to metadata/srr_list_clean.txt.
SRR_LIST="$PROJECT_DIR/metadata/srr_list_clean.txt"

FASTQ_DIR="$PROJECT_DIR/results/fastq"
FASTQC_DIR="$PROJECT_DIR/results/fastqc"
MULTIQC_DIR="$PROJECT_DIR/results/multiqc"
SALMON_DIR="$PROJECT_DIR/results/salmon"
INDEX="$PROJECT_DIR/reference/salmon_index/grch38_index"

THREADS_FASTQ=4
THREADS_SALMON=8

echo "Starting RNA-seq batch pipeline..."
echo "Project: $PROJECT_DIR"
echo "SRR list: $SRR_LIST"
echo "FASTQ dir: $FASTQ_DIR"
echo "Salmon index: $INDEX"
echo "======================================"

mkdir -p "$FASTQ_DIR" "$FASTQC_DIR" "$MULTIQC_DIR" "$SALMON_DIR" "$PROJECT_DIR/logs"

if [[ ! -f "$SRR_LIST" ]]; then
    echo "ERROR: SRR list not found: $SRR_LIST"
    exit 1
fi

if [[ ! -d "$INDEX" ]]; then
    echo "ERROR: Salmon index not found: $INDEX"
    exit 1
fi

while IFS= read -r SRR || [[ -n "$SRR" ]]; do

    SRR="${SRR//$'\r'/}"

    if [[ -z "$SRR" ]]; then
        continue
    fi

    echo "======================================"
    echo "Processing sample: $SRR"
    echo "======================================"

    FQ1_FASTQ="$FASTQ_DIR/${SRR}_1.fastq"
    FQ2_FASTQ="$FASTQ_DIR/${SRR}_2.fastq"
    FQ1_GZ="$FASTQ_DIR/${SRR}_1.fastq.gz"
    FQ2_GZ="$FASTQ_DIR/${SRR}_2.fastq.gz"

    # Decide whether to use uncompressed or compressed FASTQ
    if [[ -f "$FQ1_GZ" && -f "$FQ2_GZ" ]]; then
        FQ1="$FQ1_GZ"
        FQ2="$FQ2_GZ"
        echo "Compressed FASTQ already exists for $SRR"
    elif [[ -f "$FQ1_FASTQ" && -f "$FQ2_FASTQ" ]]; then
        FQ1="$FQ1_FASTQ"
        FQ2="$FQ2_FASTQ"
        echo "Uncompressed FASTQ already exists for $SRR"
    else
        echo "FASTQ missing for $SRR"
        echo "Downloading SRA with prefetch..."
        cd "$PROJECT_DIR"
        prefetch "$SRR"

        echo "Converting SRA to paired FASTQ..."
        fasterq-dump "$SRR" \
            --split-files \
            --threads "$THREADS_FASTQ" \
            --outdir "$FASTQ_DIR"

        FQ1="$FQ1_FASTQ"
        FQ2="$FQ2_FASTQ"
    fi

    # FastQC
    if [[ ! -f "$FASTQC_DIR/${SRR}_1_fastqc.zip" || ! -f "$FASTQC_DIR/${SRR}_2_fastqc.zip" ]]; then
        echo "Running FastQC for $SRR..."
        fastqc "$FQ1" "$FQ2" -o "$FASTQC_DIR"
    else
        echo "FastQC already exists for $SRR"
    fi

    # Salmon quantification
    if [[ -f "$SALMON_DIR/$SRR/quant.sf" ]]; then
        echo "Salmon quant already exists for $SRR. Skipping."
    else
        echo "Running Salmon quantification for $SRR..."
        salmon quant \
            -i "$INDEX" \
            -l A \
            -1 "$FQ1" \
            -2 "$FQ2" \
            -p "$THREADS_SALMON" \
            -o "$SALMON_DIR/$SRR"
    fi

    # Compress FASTQs after successful QC + Salmon
    if [[ -f "$FQ1_FASTQ" && -f "$FQ2_FASTQ" ]]; then
        echo "Compressing FASTQ files for $SRR to save disk space..."
        gzip "$FQ1_FASTQ" "$FQ2_FASTQ"
    else
        echo "FASTQ files already compressed for $SRR"
    fi

    # Remove SRA archive folders after successful processing
    echo "Removing temporary SRA archive folders for $SRR if present..."
    rm -rf "$PROJECT_DIR/$SRR" "$FASTQ_DIR/$SRR" "$PROJECT_DIR/data/raw_fastq/$SRR"

    echo "DONE: $SRR"

done < "$SRR_LIST"

echo "======================================"
echo "Running MultiQC summary..."
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --force

echo "ALL SAMPLES COMPLETED SUCCESSFULLY"

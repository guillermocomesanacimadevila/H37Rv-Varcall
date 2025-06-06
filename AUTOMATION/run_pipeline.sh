#!/bin/bash

set -euo pipefail

echo "==========================================="
echo "     🧬 Variant Calling Pipeline Runner      "
echo "==========================================="

read -rp "📁 Path to FASTA reference file: " REF
read -rp "📁 Path to FASTQ read 1: " READ1
read -rp "📁 Path to FASTQ read 2 (leave empty if SE): " READ2
read -rp "📝 Output prefix (e.g. sample1): " OUTPREFIX

nextflow run main.nf \
  --ref "$REF" \
  --reads1 "$READ1" \
  --reads2 "$READ2" \
  --outprefix "$OUTPREFIX"

echo "✅ Pipeline complete. See Results/ folder."

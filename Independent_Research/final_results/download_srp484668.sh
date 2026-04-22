#!/bin/bash
# download_srp484668.sh
# Downloads all paired-end FASTQs for Kim & Lee 2024 (SRP484668)
# Run this on the login node (needs internet access)
# Usage: bash download_srp484668.sh

set -euo pipefail

PROJ_DIR=/export/home/bio_class/barobinson4742/Bioinformatics
cd $PROJ_DIR
mkdir -p sra

# ---- 1. Fetch sample metadata from ENA ----
echo "Fetching sample list from ENA..."
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP484668&result=read_run&fields=run_accession,library_layout,read_count,base_count,sample_title,fastq_ftp&format=tsv" \
    > srp484668_metadata.tsv

# Extract just the run accessions (paired-end only) into a samples file
awk -F'\t' 'NR>1 && $2=="PAIRED" {print $1}' srp484668_metadata.tsv > samples.txt

TOTAL=$(wc -l < samples.txt)
echo "Found $TOTAL paired-end samples"

# ---- 2. Download FASTQs from ENA via FTP ----
# Uses the direct FASTQ URLs from ENA (no SRA toolkit needed)
# Downloads are backgrounded in batches to avoid overwhelming the connection

COUNT=0
while IFS=$'\t' read -r SRR LAYOUT READS BASES TITLE FTP_URLS; do
    # Skip header
    [ "$SRR" = "run_accession" ] && continue
    # Skip non-paired
    [ "$LAYOUT" != "PAIRED" ] && continue

    COUNT=$((COUNT + 1))

    # Check if already downloaded
    if [ -f "sra/${SRR}_1.fastq.gz" ] && [ -f "sra/${SRR}_2.fastq.gz" ]; then
        echo "[$COUNT/$TOTAL] $SRR already downloaded, skipping"
        continue
    fi

    echo "[$COUNT/$TOTAL] Downloading $SRR ($TITLE)..."

    # FTP_URLS is semicolon-separated: ftp.sra.ebi.ac.uk/vol1/fastq/...
    URL1=$(echo "$FTP_URLS" | tr ';' '\n' | grep '_1.fastq.gz' | head -1)
    URL2=$(echo "$FTP_URLS" | tr ';' '\n' | grep '_2.fastq.gz' | head -1)

    if [ -n "$URL1" ] && [ -n "$URL2" ]; then
        wget -q -O "sra/${SRR}_1.fastq.gz" "https://$URL1" &
        wget -q -O "sra/${SRR}_2.fastq.gz" "https://$URL2" &

        # Limit to 4 concurrent downloads to be a good network citizen
        if (( COUNT % 4 == 0 )); then
            wait
            echo "  Batch complete, continuing..."
        fi
    else
        echo "  WARNING: Could not find FTP URLs for $SRR, skipping"
    fi

done < srp484668_metadata.tsv

# Wait for any remaining downloads
wait
echo "===== All downloads complete ====="
echo "Verify with: ls sra/*.fastq.gz | wc -l  (should be $((TOTAL * 2)))"

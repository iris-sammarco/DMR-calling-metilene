#!/bin/bash

#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=40gb:scratch_local=15gb
#PBS -N prepare_metilene_input
#PBS -j oe

## Author: Iris Sammarco
## Date: 04/2025
## Purpose:
# Convert Bismark coverage files filtered for coverage >= 5 to bed files with absolute methylation ratios,
# then generate merged unionbed files per methylation context and treatment comparison for metilene DMR (differentially methylated regions) calling.
#
# Inputs:
# - Directory containing Bismark filtered coverage files
# - These files have columns: chromosome, start, end, methylation %, count methylated, count unmethylated
# - These files are expected to be named as: <context>_<sample_suffix>.bed, for example: CpG_LC1.bed, CHH_SD3.bed
#
# Outputs:
# - Bed files: methylation ratio files with 4 columns [chrom, start, end, ratio] in a tmp folder
# - Metilene input unionbed files per context and comparison
#
# Requirements:
# - bedtools installed and in PATH
# - conda environment with mamba and tools activated (modify as needed)
#
# Usage:
# qsub 2_prepare_metilene_input.sh

# --- User customization ---

#COV_DIR="/path/to/bedgraph_coverage_5"
#OUTPUT_DIR="/path/to/dmrs/metilene_input"

#####
export TMPDIR=$SCRATCHDIR
# Load mamba
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate tools
COV_DIR="/storage/pruhonice1-ibot/home/irissammarco/trash"
OUTPUT_DIR="/storage/pruhonice1-ibot/home/irissammarco/trash/metilene_input"
#####

mkdir -p "${OUTPUT_DIR}/tmp_bed"

# Define methylation contexts
contexts=("CpG" "CHG" "CHH") # customize these names

# Specify methylation contexts to use (e.g., CpG, CHG, CHH)
contexts=("CpG" "CHG" "CHH")

# Define sample suffixes exactly as they appear in bed filenames (after context_ prefix)
group_LC=("LC1" "LC2" "LC3")
group_LD=("LD1" "LD2" "LD3")
group_SC=("SC1" "SC2" "SC3")
group_SD=("SD1" "SD2" "SD3")

# Define comparisons as pairs of group variable names
declare -A comparisons
comparisons["LC_LD"]="group_LC group_LD"
comparisons["SC_SD"]="group_SC group_SD"
comparisons["LC_SC"]="group_LC group_SC"
comparisons["LD_SD"]="group_LD group_SD"

# --- Step 1: Convert coverage to bed with methylation ratio and 0-based start coordinate ---
for cov_file in ${COV_DIR}/*.bed; do # if the filenames' extension is different from .bed, change it here and in the line below
    base=$(basename "$cov_file" .bed)
    awk '($5 + $6) > 0 {
        ratio = $5 / ($5 + $6);
        if ($2 == $3) $2 = $2 - 1;     # Convert to 0-based start coordinate (BED format requirement)
        printf "%s\t%s\t%s\t%.4f\n", $1, $2, $3, ratio;
    }' OFS='\t' "$cov_file" > "${OUTPUT_DIR}/tmp_bed/${base}.bed"
done

# --- Step 2: Create metilene unionbed input files for all context/comparison combinations ---
for context in "${contexts[@]}"; do
  for comparison_key in "${!comparisons[@]}"; do
    
    # Extract the two group array names from the comparisons associative array
    read -r group1_name group2_name <<< "${comparisons[$comparison_key]}"
    
    # Bash 'nameref' to access arrays by their names stored as strings in comparison map
    declare -n group1=${group1_name}
    declare -n group2=${group2_name}

    # Build input file lists by appending sample bed file paths dynamically
    in1_files=""
    for sample in "${group1[@]}"; do
      in1_files+=" ${OUTPUT_DIR}/tmp_bed/${context}_${sample}.bed"
    done

    in2_files=""
    for sample in "${group2[@]}"; do
      in2_files+=" ${OUTPUT_DIR}/tmp_bed/${context}_${sample}.bed"
    done

    output_file="${OUTPUT_DIR}/${context}_${comparison_key}_union.bed"

    # Run unionbedg with missing values filled by '-' 
    bedtools unionbedg -filler - -i $in1_files $in2_files |
    cut -f1,3- |    # Remove original start coordinate as BED start (now field 2 is union start)
    awk '{
        printf "%s\t%s", $1, $2;
        for(i=3;i<=NF;i++) {
            if($i=="-") {
              printf "\t%s", $i
            } else {
              printf "\t%.2f", $i
            }
        }
        printf "\n"
    }' |
    sort -k1,1 -k2,2n -T "${OUTPUT_DIR}/tmp_bed" > "$output_file"

  done
done

# --- Step 3: Add column header required by Metilene and format output files ---
cd "$OUTPUT_DIR" || exit

for file in *_union.bed; do
    base=$(basename "$file" _union.bed)
    # Extract context and comparison (two group names) from filename
    context=${base%%_*}
    groups=${base#*_}
    IFS="_" read -r group1 group2 <<< "$groups"

    # Construct header with sample replicate names per group (assumes 3 replicates)
    header="chrom\tpos"
    for i in {1..3}; do
        header="${header}\t${group1}_${i}"
    done
    for i in {1..3}; do
        header="${header}\t${group2}_${i}"
    done

    tmp="${file}.tmp"
    mv "$file" "$tmp"
    echo -e "$header" > "$file"
    cat "$tmp" >> "$file"
    rm "$tmp"
done
# DMR-calling-metilene
A bioinformatics pipeline for Differentially Methylated Region (DMR) calling using Metilene. Includes scripts for filtering bedgraph files, preparing Metilene input data, running Metilene, and filtering output results.

# Pipeline Overview
The pipeline consists of three main scripts executed in order:

1_filter_bedgraphs_cov5.bash

Filters input bedgraph files to retain sites with coverage above a defined threshold (default 5 reads). Produces filtered bedgraph files as output.

2_prepare_metilene_input.bash

Converts filtered bedgraph files into Metilene-compatible input files. This includes formatting and merging methylation data across samples.

3_run_metilene_and_filter.bash

Runs Metilene DMR caller on the prepared input files. Filters Metilene output to retain high-confidence DMRs based on predefined criteria.

# Input and Output
Input files: Raw bedgraph files for each sample containing methylation and coverage information.

Intermediate files: Filtered bedgraphs and formatted Metilene input files.

Output files: Filtered DMR lists called by Metilene.

Adjust parameters inside each script as needed for coverage thresholds, input file patterns, or filtering criteria.

#!/bin/bash

# # Salmon index # #
# source activate salmon1
# transcriptome="/home/ls/rachelcw/projects/salmon/transcriptome.fa"
# ~/miniconda3/envs/salmon1/bin/salmon index -t "$transcriptome" -i ds_index -k 1 
# conda deactivate

# List of BAM files
# bam_files=("/data01/private/data/cllmap/data/rna/bams/DKFZ-CLL65.Aligned.sortedByCoord.out.bam")  # Add your BAM file paths here

# Output directory for FASTQ and Salmon quant results
output_dir="quant_results"
mkdir -p "$output_dir"

# Iterate through each BAM file
# for file in `ls /private/data/cllmap/data/rna/bams/`; do
# Convert BAM to paired-end FASTQ using Picard
file="/data01/private/data/cllmap/data/rna/bams/DKFZ-CLL65.Aligned.sortedByCoord.out.bam"
base_name=$(basename "$file" .bam)
fastq_output="$output_dir/${base_name}.fastq"
echo "Processing $file"
~/miniconda3/envs/picard/bin/picard SamToFastq I="$file" FASTQ="$fastq_output"
echo "Finished converting $file to FASTQ"
# Run Salmon quantification on the FASTQ
salmon_output="$output_dir/${base_name}_quant"
~/miniconda3/envs/salmon1/bin/salmon quant -i ds_index -l A -r "$fastq_output" -p 8 -o "$salmon_output"
echo "Finished quantifying $file"
# Remove the FASTQ file
rm "$fastq_output"
echo "Finished processing $file"
# done

# # Aggregate Salmon quant results
# # aggregate_output="$output_dir/aggregate_quant"
# # salmon quantmerge --quants "$output_dir"/*_quant --output "$aggregate_output"

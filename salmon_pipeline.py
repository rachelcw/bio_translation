import os
import subprocess
import glob
from concurrent.futures import ProcessPoolExecutor as Pool

# GLOBAL VARIABLES
input_dir = "/private/data/cllmap/data/rna/bams/"
# Specify Salmon index and output directory
salmon_index = "/home/ls/rachelcw/projects/salmon/ds_fasta_index"
output_dir = "quant_results"
def convert_bam_fastq(bam_path, fastq_prefix):
    subprocess.run(["source activate picard"], shell=True)
    subprocess.run([
        "/home/ls/rachelcw/miniconda3/envs/picard/bin/picard",
        "SamToFastq",
        "INPUT=" + bam_path,
        "FASTQ=" + fastq_prefix + "_R1.fastq",
        "SECOND_END_FASTQ=" + fastq_prefix + "_R2.fastq",
    ], check=True)
    print("SamToFastq completed.")

# Loop through each BAM file
# for bam_file in bam_files:
def pipeline(bam_file):
    bam_path = os.path.join(str(input_dir), str(bam_file))
    # Convert BAM to paired-end FASTQ using Picard's SamToFastq
    output_prefix = os.path.splitext(os.path.basename(bam_file))[0]
    fastq_prefix = os.path.join(output_dir, output_prefix)
    print("Processing " + output_prefix)
    # convert_bam_fastq(bam_path, fastq_prefix)
    # fastq_done=[os.path.splitext(os.path.basename(bam_file))[0] for f in os.listdir(output_dir) if f.endswith(".fastq")]
    
    # Run Salmon quantification on the generated FASTQ files
    salmon_output_dir = os.path.join(output_dir, output_prefix + "_quant")
    subprocess.run([
        "/home/ls/rachelcw/miniconda3/envs/salmon1/bin/salmon",
        "quant",
        "-i", salmon_index,
        "-l", "A",
        "-1", fastq_prefix + "_R1.fastq",
        "-2", fastq_prefix + "_R2.fastq",
        "-p", "80",  # Number of threads
        "-o", salmon_output_dir
    ], check=True)
    print("Salmon quantification completed.")
    # Remove the generated FASTQ files
    os.remove(fastq_prefix + "_R1.fastq")
    os.remove(fastq_prefix + "_R2.fastq")
    print("fastq files removed")
    print("Finished processing " + output_prefix)


if __name__== "__main__":
    bam_files = [f for f in os.listdir(input_dir) if f.endswith(".bam")]
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # Run with multiprocessing
    with Pool(max_workers=50) as pool:
        pool.map(pipeline, bam_files)

    # Aggregate Salmon quant results
    aggregate_output_dir = os.path.join(output_dir, "aggregate_results")
    os.makedirs(aggregate_output_dir, exist_ok=True)

    # Loop through the quant result directories and copy the quant.sf files
    for subdir in os.listdir(output_dir):
        if subdir.endswith("_quant"):
            quant_sf = os.path.join(output_dir, subdir, "quant.sf")
            new_quant_sf = os.path.join(aggregate_output_dir, subdir + "_quant.sf")
            os.rename(quant_sf, new_quant_sf)

    print("Quantification and aggregation completed.")

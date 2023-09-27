import os
import subprocess
import glob
from concurrent.futures import ProcessPoolExecutor as Pool

# # Salmon index # #
# source activate salmon1
# transcriptome="/home/ls/rachelcw/projects/salmon/transcriptome.fa"
# ~/miniconda3/envs/salmon1/bin/salmon index -t "$transcriptome" -i ds_index


# GLOBAL VARIABLES
# input_dir = "/private/data/cllmap/data/rna/bams/"
input_dir = "/data01/private/data/cllmap/rnaseq/bams_cll11/"

# Specify Salmon index and output directory
salmon_index = "/home/ls/rachelcw/projects/salmon/ds_fasta_index"
output_dir = "quant_results"

# Convert BAM to paired-end FASTQ using Picard's SamToFastq
# check that fastq file has 4 lines per read
# check that you have enoungh space in JVM, if not, 
    # increase the memory in picard file in /home/ls/rachelcw/miniconda3/envs/picard/bin/picard in the part of java -Xmx
def convert_bam_fastq(bam_path, fastq_prefix):
    # subprocess.run(["source activate picard"], shell=True)
    subprocess.run([
        "/home/ls/rachelcw/miniconda3/envs/picard/bin/picard",
        "SamToFastq",
        "INPUT=" + bam_path,
        "FASTQ=" + fastq_prefix + "_R1.fastq",
        "SECOND_END_FASTQ=" + fastq_prefix + "_R2.fastq",
    ], check=True)
    print("SamToFastq completed.")

# Run Salmon quantification on the generated FASTQ files
# check that you
def sam_quant(fastq_prefix):
    salmon_output_dir = fastq_prefix + "_quant"
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

# when you have a list of fastq files already and you want to run salmon on them
def fastq_to_sam():
    fastq_done=[os.path.splitext(os.path.basename(f))[0] for f in os.listdir(output_dir) if f.endswith(".fastq")]
    fastq_done=[fastq.replace(r'_R1','') for fastq in fastq_done ]
    fastq_done=[fastq.replace(r'_R2','') for fastq in fastq_done ]
    fastq_done=set(fastq_done)
    for fq_prefix in fastq_done:
    # Run Salmon quantification on the generated FASTQ files
        print("Processing " + fq_prefix)
        salmon_output_dir = os.path.join(fq_prefix + "_quant")
        print(salmon_output_dir)
        fastq_prefix = os.path.join(output_dir, fq_prefix)
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

# Loop through each BAM file
def pipeline(bam_file):
    bam_path = os.path.join(str(input_dir), str(bam_file)) 
    print(bam_path)
    # # Convert BAM to paired-end FASTQ using Picard's SamToFastq
    output_prefix = os.path.splitext(os.path.basename(bam_file))[0]
    fastq_prefix = os.path.join(output_dir, output_prefix)

    print("Processing " + fastq_prefix)
    convert_bam_fastq(bam_path, fastq_prefix)
    sam_quant(fastq_prefix)
    # # # Remove the generated FASTQ files
    os.remove(fastq_prefix + "_R1.fastq")
    os.remove(fastq_prefix + "_R2.fastq")
    print("fastq files removed")
    print("Finished processing " + fastq_prefix)

def gcll_proteomics_bam_not_quant_list():
    # Open the text file for reading
    file_path = "/home/ls/rachelcw/projects/salmon/gcll_proteomics_bam_not_quant.txt" 
    with open(file_path, "r") as file:
        content = file.read()

    # Split the content using commas and create a list
    elements_list = content.split(",")

    # Remove leading and trailing whitespaces from each element
    elements_list = [element.strip() for element in elements_list]


    return elements_list


if __name__== "__main__":
     
    # Create the output directory if it doesn't exist
    # os.makedirs(output_dir, exist_ok=True)

    # # Loop through each BAM file and run the pipeline 
    bam_files = [f for f in os.listdir(input_dir) if f.endswith(".bam")]
    # prefix_bam_files=[os.path.splitext(os.path.basename(f))[0].replace('.out','') for f in bam_files]
    
    bams_not_quant=gcll_proteomics_bam_not_quant_list() # startwith
    list_samples = [element for element in bam_files if any(element.startswith(prefix) for prefix in bams_not_quant)]
    print(len(list_samples))
    # quant_files=[f for f in os.listdir(output_dir) if f.endswith("_quant")]
    # sf_files=[f for f in quant_files for sf in os.listdir(os.path.join(output_dir,f)) if sf.endswith(".sf") ]
    # prefix_sf_files=[os.path.splitext(os.path.basename(f))[0] for f in sf_files]
    # list_samples=[f'{s}.out.bam' for s in prefix_bam_files if s not in prefix_sf_files ]
    # print(list_samples)

    # # run single process
    # pipeline(list_samples[0])

    # # run multi process
    with Pool(max_workers=80) as pool:
        pool.map(pipeline, list_samples)

    

    # # Aggregate Salmon quant results
    # aggregate_output_dir = os.path.join(output_dir, "aggregate_results")
    # os.makedirs(aggregate_output_dir, exist_ok=True)

    # # Loop through the quant result directories and copy the quant.sf files
    # for subdir in os.listdir(output_dir):
    #     if subdir.endswith("_quant"):
    #         quant_sf = os.path.join(output_dir, subdir, "quant.sf")
    #         new_quant_sf = os.path.join(aggregate_output_dir, subdir + "_quant.sf")
    #         os.rename(quant_sf, new_quant_sf)

    # print("Quantification and aggregation completed.")

import sys
import numpy as np
from optparse import OptionParser
import pyfaidx as fa
import pyensembl
import csv
from Bio.Seq import reverse_complement
from Bio.Seq import translate



# codon table to prtoein #
table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    
def convert_exon_pos_to_dict(exons_file):
    """
    the function get a file with the start and end positions of exons and return dict of chr and its exons' positions
    
    
    Args:
        exons_file (txt file)

    Returns:
        list of the position of each exon[[start,end]...]
    """
    exons=dict()
    final_position=dict()
    with open(exons_file) as f:
        # Read the file as a CSV
        reader = csv.reader(f, delimiter="\t") #TODO what is sep
        # Iterate over the rows in the file
        for row in reader:
            # Get the start and end positions
            chr=row[0]
            if chr.startswith('chr'):
               chr=chr.replace('chr',"")
            start = int(row[1])
            end = int(row[2])
            if chr not in exons.keys():
                exons[chr]=[[start,end]]
            else:
                exons[chr].append([start,end])
        for chr in exons.keys():
            exons[chr].sort()
            final_position[chr]=[np.min(exons[chr]),np.max(exons[chr])]
    return final_position,exons


def get_seq(exons,fasta):
    exons_seq=dict()
    """
    Extract the exons sequence from fasta file 
    Args:
        exons (list): the list from the function convert_exon_pos_to_list
        fasta (fasta file): fasta file of genes

    Returns:
        list with the sequences of the exons for each gene from fasta file 
    """
    for chr in exons.keys():        
        sequence=""
        for start,end in exons[chr]:
            seq=fasta[chr][start:end].seq
            sequence=sequence+seq
        exons_seq[chr]=sequence
    return exons_seq


def start_read(rna_seq,n,start):
    """
    Edit rna seq to start position n
    Change the rna seq to start from AUG/5/3 prime
    Args:
        rna_seq (string): rna seq of gene 
        n (int): the start position to read the seq
        start (int): how to read the seq
    Returns:
        seq: the sequence after editing
    """
    seq=rna_seq[n:] # TODO it's work the same for all sequnces- change?
    if start == 0:
        #start from ATG
            seq=seq[seq.find('ATG'):]
    elif start == 3:
        #start from 3 prime
            seq=seq[::-1] #reverse
    elif start == 5:
            seq=seq
    return seq

def frame_read(seq,frame):
    """
    The function get sequence and the reading frame and return the translated seq
    Args:
        seq (string): the sequence for translating
        frame (int): the reading frame

    Returns:
        list of translated seq to protein by the reading frame 
    """
    translate=[]
    if frame == 1:
            translate.append(translation(seq))
    if frame == 3 or frame == 6:
           frame1=translation(seq)
           frame2=translation(seq[1:])
           frame3=translation(seq[2:])
           translate.extend([frame1,frame2,frame3])
    if frame == 6:
            frame4=translation(seq[::-1])
            frame5=translation(seq[::-2])
            frame6=translation(seq[::-3])
            translate.extend([frame4,frame5,frame6])
    return translate


def translation(seq):
    """
    The function get rnaseq and translates to protein by the codon table
    Args:
        seq (string): the sequence

    Returns:
        string of amino acid
    """
    seq=transcribe(seq) # transcribe to cdna(rna) for translation
    protein = []
    for i in range(0,len(seq),3):
        codon = seq[i:i+3]
        if codon in table:
            aminoacid = table[codon]
            # print(aminoacid)
            if aminoacid== '*':
                 # stop translation when there is a stop codon 
                break
            protein.append(aminoacid)
        else:
            protein.append("N")         
    return "".join(protein)
    # protein=translate(seq)
    # return protein

# def transcript_gtf(fasta):
#     genes=[]
#     # Open the GTF file
#     with open("/private1/private/resources/gencode19_noChrPrefix_mitoMT.gtf") as f:
#         # Read the GTF file as a CSV
#         reader = csv.reader(f, delimiter="\t")

#         # Iterate over the rows in the GTF file
#         for row in reader:
#             if row[0].startswith('#'):
#                 continue
#             # Check if the feature is a transcript
#             elif row[2] == "transcript":
#                 # Get the chromosome and start/end positions of the transcript
#                 chromosome = row[0]
#                 start = int(row[3]) - 1  # GTF is 1-based, pyfaidx is 0-based
#                 end = int(row[4])

#                 # Get the sequence object for the chromosome
#                 seq = fasta[chromosome]

#                 # Extract the transcript sequence
#                 transcript_seq = seq[start:end]
#                 genes.append([chromosome,start,end,transcript_seq])

# def exon_gtf():
# # Open the GTF file
#     with open("path/to/file.gtf") as f:
#     # Read the GTF file as a CSV
#         reader = csv.reader(f, delimiter="\t")

#         # Iterate over the rows in the GTF file
#         for row in reader:
#             # Check if the feature is an exon
#             if row[2] == "exon":
#                 # Get the chromosome and start/end positions of the exon
#                 chromosome = row[0]
#                 start = int(row[3]) - 1  # GTF is 1-based, pyfaidx is 0-based
#                 end = int(row[4])

#                 # Get the sequence object for the chromosome
#                 seq = fasta[chromosome]

#                 # Extract the exon sequence
#                 exon_seq = seq[start:end]

#                 # Print the exon sequence
#                 print(exon_seq)

def transcribe(seq):
    # Convert the transcript sequence to cDNA
    return reverse_complement(seq)

def information_gtf(chr,start,end):
    data = pyensembl.Genome(
    reference_name='GRCh37',
    annotation_name='my_genome_lab',
    gtf_path_or_url='/private1/private/resources/gencode19_noChrPrefix_mitoMT.gtf')
    data.index()
    gene_name=data.gene_names_at_locus(chr,int(start),int(end))
    if gene_name==[]:
        gene_name='intergenetic'
    return gene_name

if __name__== "__main__":
    parser = OptionParser()
    parser.add_option("-i","--exons",dest="input",
                help="txt file with start and end positions of the exons")
    parser.add_option("-a","--fasta", dest="fasta",
                help="fasta file- *.fa")
    # parser.add_option("ifasta",
    #             help="fasta index file- *.fai")
    parser.add_option("-n", "--start-position", dest="n",default='0',
                help="0- start from start postion(default), n- start read the sequence from the n position ") #TODO
    parser.add_option("-s", "--start", dest="start",choices=['0','5','3'], default='0',
                help="0- start from start codon(default) 5- start from 5 prime 3-start from 3 prime ")
    parser.add_option("-f", "--frame", dest="frame",default='1',choices=['1','3','6'],
                help="read frame: 1(default), 3 (1,2,3 frame), 6 (1,2,3 frame from 5 and 3 prime")
    parser.add_option("-o", "--out", dest="output", default = 'protein',
                  help="output filename with its path ")
    
    (options, args) = parser.parse_args()
    
    if options.input == None:
        sys.stderr.write("Error: no input file provided...\n")
        exit(0)
        
    if options.fasta == None:
        options.fasta="/private1/private/resources/Homo_sapiens_assembly19.fasta"
        
    # if options.ifasta == None:
    #     options.ifasta="/private1/private/resources/Homo_sapiens_assembly19.fasta.fai"

    exons_file=options.input
    fasta =fa.Fasta(options.fasta)
    # fai=options.ifasta
    n=int(options.n)
    startr=int(options.start)
    frame=int(options.frame)
    final_position,exons=convert_exon_pos_to_dict(exons_file) #dict[chr]=[[start,end]...]
    genes_seq=get_seq(exons,fasta) #dict[chr]=seq
    with open(options.output, "w") as f:
        for chr,seq in genes_seq.items():
            seq=start_read(seq,n, startr)
            translate=frame_read(seq, frame)
            # strand,gene_name=information_gtf(chr,final_position[chr][0],final_position[chr][1])
            gene_name=information_gtf(chr,final_position[chr][0],final_position[chr][1])
            title=f'>{chr}:{final_position[chr][0]}-{final_position[chr][1]}|{gene_name}|'
            for i,protein in enumerate(translate):
                if i <=2:
                    f.write(title)
                    f.write(f'frame +{i+1}\n')
                elif i>2:
                    f.write(title)
                    f.write(f'frame -{i-2}\n')
                f.write(protein+'\n')
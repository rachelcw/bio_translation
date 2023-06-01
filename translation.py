import sys
import numpy as np
from optparse import OptionParser
import pyfaidx as fa
import pyensembl
import csv
import pandas as pd
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

def get_args():
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
    
    options, args= parser.parse_args()

    if options.input == None:
        sys.stderr.write("Error: no input file provided...\n")
        exit(0)
        
    if options.fasta == None:
        options.fasta="/private1/private/resources/Homo_sapiens_assembly19.fasta"

    # parser input argument 
    input_file=options.input
    fasta =fa.Fasta(options.fasta)
    n=int(options.n)
    startr=int(options.start) # strat read from..
    frame=int(options.frame)
    if frame>1:
        # when the frame are not one, start read from the beginning, not from AUG
        startr=5
    
    return input_file,fasta,n,startr,frame,options.output
    
    
# def get_junction_information(input_file):
#     """
#     the function get a file with the start and end positions of exons and return dict of chr and its exons' positions
    
#     Args:
#         exons_file (txt file)

#     Returns:
#         list of the position of each exon[[start,end]...]
#     """
#     # bed_file = collections.namedtuple('junction', ['name', 'age', 'DOB'])
#     # transcript=collections.namedtuple("transcript",['chr', 'start', 'end','strand','gene_id'])
#     junction_df=pd.read_csv(input_file,sep='\t',header=None,names=["chr","start","end","strand","ENST","ENSG","junction"])
#     junc_tranacript=dict()
#     for junction in junction_df["junction"].unique():
#         junc_tranacript[junction]=[t for t in junction_df.loc[junction_df["junction"]==junction,"ENST"].unique()]

#     for junction,transcripts in junc_tranacript.items():
#         for t_id in transcript:
            
#         with open(input_file) as f:
#         # Read the file as a CSV
#         reader = csv.reader(f, delimiter="\t") 
#         # Iterate over the rows in the file
#         for row in reader:
#             # Get the start and end positions
#             chr=row[0]
#             if chr.startswith('chr'):
#                chr=chr.replace('chr',"")
#             start = int(row[1])
#             end = int(row[2])
#             strand=row[3]
#             t_id=row[4]
#             gene_id=row[5]
#             junction=row[6]
#             transcript=transcript(chr,start,end,strand,gene_id)
#             if junction not in junction_dict.keys():
#                 junction_dict[junction]={t_id:transcript(chr,start,end,strand,gene_id)}
                
#             else:
#                 if t_id not in junction_dict[junction].keys():
#                     junction_dict[junction]={t_id:transcript(chr,start,end,strand,gene_id)}
#                 else:
#                     junction_dict[junction][t_id].append(transcript(chr,start,end,strand,gene_id))
#     return junction_dict


def get_seq(locus,fasta):
    # junction_seq=dict()
    """
    Extract the exons sequence from fasta file 
    Args:
        exons (list): the list from the function convert_exon_pos_to_list
        fasta (fasta file): fasta file of genes

    Returns:
        list with the sequences of the exons for each gene from fasta file 
    """
    sequence=""
    for exon in locus:        
        seq=fasta[exon.chr][exon.start:exon.end].seq
        sequence=sequence+seq
    return sequence
    # for junc in junction_dict.keys():        
    #     sequence=""
    #     for transrcript in junc.keys():
    #         seq=fasta[exon.chr][exon.start:exon.end].seq
    #         sequence=sequence+seq
    #     junction_seq[junc]=sequence
    # return junction_seq


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
           print("frame1", frame1)
           frame2=translation(seq[1:])
           print("frame2", frame2)
           frame3=translation(seq[2:])
           print("frame3", frame3)
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
        elif len(codon)%3 == 0:
            protein.append("N")         
    return "".join(protein)
    

def information_gtf(chr,start,end):
    data = pyensembl.Genome(
    reference_name='GRCh37',
    annotation_name='my_genome_lab',
    gtf_path_or_url='/private1/private/resources/gencode19_noChrPrefix_mitoMT.gtf')
    data.index()
    gene_name=data.gene_names_at_locus(chr,int(start),int(end))
    if gene_name==[]:
        gene_name='intergenic'
    return gene_name
    
# def create_output_file(junction_dict,genes_seq,n,startr,frame,output_path):
#      with open(output_path, "w") as f:
#         for junc in junction_dict.keys():
#             seq=genes_seq[junc]
#         for chr,seq in genes_seq.items():
#             seq=start_read(seq,n, startr)
#             translate=frame_read(seq, frame)
#             # strand,gene_name=information_gtf(chr,final_position[chr][0],final_position[chr][1])
#             gene_name=information_gtf(chr,final_position[chr][0],final_position[chr][1])
#             title=f'>{chr}:{final_position[chr][0]}-{final_position[chr][1]}|{gene_name}|'
#             for i,protein in enumerate(translate):
#                 if i <=2:
#                     f.write(title)
#                     f.write(f'frame+{i+1}\n')
#                 elif i>2:
#                     f.write(title)
#                     f.write(f'frame-{i-2}\n')
#                 f.write(protein+'\n')

def get_intron(junc):
      import re
      return re.sub(r'clu_\d+_', '', junc)

def write_output(seq,gene_id,t_id,junc,file):
    with open(file, "w") as f:
        f.write(f'>{gene_id}|{t_id}|{junc}\n')
        f.write(seq+'\n')

if __name__== "__main__":
    input_file,fasta,n,startr,frame,output_path=get_args()
    junction_df=pd.read_csv(input_file,sep='\t',header=None,names=["chr","start","end","strand","ENST","ENSG","junction"])
    with open(output_path, "w") as f:
        for junc in junction_df["junction"].unique():
            transcript=junction_df.loc[junction_df["junction"]==junc].unique()
            for t in transcript:
                exons=transcript.loc[transcript["ENST"]==t]
                locus=[(exon.chr,exon.start,exon.end) for exon in exons.itertuples()]
                seq=get_seq(locus,fasta)
                intron=get_intron(junc)
                f.write(f'>{exons["ENSG"][0]}|{t}|{intron}\n')
                f.write(seq+'\n')
             
        

    # junction_dict=get_junction_information(input_file) #dict[chr]=[[start,end]...]
    # genes_seq=get_seq(junction_dict,fasta) #dict[chr]=seq #TODO what happen when we have more than one gene from the same chr
    # create_output_file(junction_dict,genes_seq,n,startr,frame,output_path)
    
    
    
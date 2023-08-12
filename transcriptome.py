# create transcriptome file
#python translation.py -i "/data01/private/projects/splicing_cll/results/proteomics/analysis.20230724/mutated_cll_sf3b1_proteomics.20230724.novel.txt" -o /home/ls/rachelcw/projects/BIO/transcriptome_mut.fa
import pandas as pd
import pyfaidx as fa
from translation import get_seq, get_args, start_read, get_intron

def create_transcriptome_only_transcripts():
    transcriptome_file_path="/home/ls/rachelcw/projects/salmon/transcriptome.fa"
    transcript_file=open(transcriptome_file_path,"w")
    with open("/home/ls/rachelcw/projects/salmon/transcriptome_full_head.fa", "r") as f:
        transcripts={}
        for line in f:
            line=line.strip()
            if line.startswith(">"):
                transcript=line.split("|")[1]
                if transcript in transcripts.keys():
                    transcripts[transcript]+=1
                else:
                    transcripts[transcript]=0
                transcript_file.write(f'>{transcript}_DSJ{transcripts[transcript]}\n')
            else:
                transcript_file.write(line+'\n')
    reference_file_path_mut="/data01/private/projects/splicing_cll/results/proteomics/analysis.20230713/mutated_cll_sf3b1_proteomics_reference.txt"
    reference_file_path_unmut="/data01/private/projects/splicing_cll/results/proteomics/analysis.20230713/unmutated_cll_sf3b1_proteomics_reference.txt"
    reference_file_mut=pd.read_csv(reference_file_path_mut,sep='\t',header=None,names=["chr","start","end","strand","ENST","ENSG"])
    reference_file_unmut=pd.read_csv(reference_file_path_unmut,sep='\t',header=None,names=["chr","start","end","strand","ENST","ENSG"])
    reference_file=pd.concat([reference_file_mut,reference_file_unmut])
    reference_file.drop_duplicates(keep='first',inplace=True)
    for gene in reference_file["ENSG"].unique():
        transcript=reference_file[reference_file["ENSG"]==gene]
        for t in transcript["ENST"].unique():
            exons=transcript[transcript["ENST"]==t]
            locus=[(exon.chr,exon.start,exon.end) for exon in exons.itertuples()]
            strand=exons["strand"].unique()[0]
            fasta="/private1/private/resources/Homo_sapiens_assembly19.fasta"
            fasta =fa.Fasta(fasta)
            seq=get_seq(locus,fasta,strand)
            transcript_file.write(f'>{t}\n')
            transcript_file.write(f'{seq}\n')
    transcript_file.close()

def merge_transcriptome_files():
    mut_file=open("/home/ls/rachelcw/projects/salmon/transcriptome_mut.fa","r")
    unmut_file=open("/home/ls/rachelcw/projects/salmon/transcriptome_unmut.fa","r")
    head_seq_dict={}
    for line in unmut_file:
        line=line.strip()
        if line.startswith(">"):
            head=line
        else: # line = seq
            if head not in head_seq_dict.keys():
                head_seq_dict[head]=line
    for line in mut_file:
        line=line.strip()
        if line not in head_seq_dict.keys():
            if line.startswith(">"):
                head=line
            else: # line = seq
                head_seq_dict[head]=line
    with open("/home/ls/rachelcw/projects/salmon/transcriptome_full_head.fa", "w") as f:
        for head in head_seq_dict.keys():
            f.write(head+'\n')
            f.write(head_seq_dict[head]+'\n')

    mut_file.close()
    unmut_file.close()

def create_transcriptome_file_by_analysis(input_file,fasta,n,startr,frame,output_path):
    with open(output_path, "w") as f:
        data=pd.read_csv(input_file,sep='\t',header=None,names=["chr","start","end","strand","ENST","ENSG","junction","start_read"])
        for junc in data["junction"].unique():
            transcript=data[data["junction"]==junc]
            gene_id=transcript["ENSG"].unique()[0]
            for t in transcript["ENST"].unique():
                exons=transcript[transcript["ENST"]==t]
                locus=[(exon.chr,exon.start,exon.end) for exon in exons.itertuples()]
                strand=exons["strand"].unique()[0]
                seq=get_seq(locus,fasta,strand)
                startr=exons["start_read"].unique()[0] # 5-from the start , 0-from start codon                
                seq=start_read(seq,n, startr)
                intron=get_intron(junc)
                f.write(f'>{gene_id}|{t}|{intron}|\n')
                f.write(seq+'\n')

if __name__== "__main__":
    # input_file,fasta,n,startr,frame,output_path=get_args()
    # create_transcriptome_file_by_analysis(input_file,fasta,n,startr,frame,output_path)
    merge_transcriptome_files()
    create_transcriptome_only_transcripts()

""" run salmon
cd salmon
salmon quant -t transcriptome.fa  -a /data01/private/data/cllmap/data/rna/bams/DKFZ-CLL65.Aligned.sortedByCoord.out.bam /data01/private/data/cllmap/data/rna/bams/DKFZ-CLL72.Aligned.sortedByCoord.out.bam --validateMappings -o ds_salmon_quant
"""
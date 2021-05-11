import streamlit as st
import sys
import SessionState
from Gene import *
from Protein_isoform import *
from Streamlit_community import *
from Input_flow import *
from Streamlit_Pop_ups import *
from Alignment import *
from Visualise_Alignment import *
from User_Input_Preparation import *
from Input_flow import *
from Table_Generation import *
from PIL import Image
from Statistics import *
from Comparison import *


#load human isoform library
list_of_gene_objects = Input_flow.import_data_from_github('list_of_gene_objects_19th_april.txt.gz')


#iterate through library
#for gene in list_of_gene_objects:
#    prepare_fasta_file_as_input_for_MSA_only ensembl forers
#    post fasta files
#    wait for retrieval of files
#
#run alignment and exon function for reference Id with needle-man wunsch
#extract alignment of fasta file batch of collection of isoforms
#run exon function for MSA aligment
#
#run comparison between the two lists and make statistics? similarity,

for gene in list_of_gene_objects:
    list_of_fastas = Comparison.extract_protein_sequences_in_fasta_from_gene(gene)
    if not list_of_fastas or len(list_of_fastas)<3:
        continue

    fasta_string = "\n".join(list_of_fastas)
    print(fasta_string)






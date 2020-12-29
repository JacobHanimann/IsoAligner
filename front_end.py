import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("TkAgg")
import matplotlib as plt
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np
from PIL import Image
import gzip
import pickle
import csv
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import base64
import re



#Streamlit website
def main():
    """ Isoform Mapping Tool """
    st.title("AminoAcid Isoform Mapper")

    def get_binary_file_downloader_html(bin_file, file_label='File', name_button='download'):
        with open(bin_file, 'rb') as f:
            data = f.read()
        bin_str = base64.b64encode(data).decode()
        href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">{name_button} {file_label}</a>'
        return href

    activity = ['Mapping tool', 'Download pre-computed data', 'About & Source Code']
    st.sidebar.markdown("### Menu")
    choice = st.sidebar.selectbox("", activity)
    st.sidebar.write("\n")
    st.sidebar.write("\n")

    if choice == 'Mapping tool':
        st.subheader("A simple tool to align isoforms globally")
        st.write("Align isoforms with the Needleman-Wunsch algorithm and set the minimal exon length to discard fasely mapped positions (random matches) of two distinct exons."
                 " The table of correctly mapped positions can be downloaded as a file in several formats. A preview of the alignments is displayed dynamically. ")
        st.write("--------------------------")
        st.markdown("#### Input")
        fasta1 = st.text_area('Paste Amino Acid sequence of reference isoform: ', '''''')
        agree = st.checkbox("upload fasta file of reference isoform?")
        if agree:
            fasta1 = st.file_uploader("Upload 1st FASTA File", type=["fasta", "fa","txt"])
        st.write("\n")
        option = st.selectbox(
            'Select alternative isoforms to align',
             ['Insert own sequence',"Available on Refseq", "Available on Ensembl", "All Available Isoforms"])
        if option == "Insert own sequence":
            fasta2 = st.text_area('Type in Amino Acid sequence of alternative isoform: ', '''''')
            agree2 = st.checkbox("upload fasta file of alternative isoform?")
            if agree2:
                fasta2 = st.file_uploader("Upload 2st FASTA File", type=["fasta", "fa","txt"])
        st.write("--------------------------")
        st.sidebar.write("\n")
        st.sidebar.markdown("#### Needleman-Wunsch-Algorithm Parameters:")
        st.sidebar.write("\n")
        match= st.sidebar.number_input("Match:", min_value=None, max_value=None, value=1, step=None, format=None, key=None)
        mismatch= st.sidebar.number_input("Mismatch:", min_value=None, max_value=None, value=-2, step=None, format=None, key=None)
        open_gap_penalty= st.sidebar.number_input("Open Gap penalty:", min_value=None, max_value=None, value=-1.75, step=None, format=None, key=None)
        gap_extension_penalty= st.sidebar.number_input("Gap Extension penalty:", min_value=None, max_value=None, value=0, step=None, format=None, key=None)
        st.sidebar.write("\n")
        st.sidebar.markdown("#### Minimal Exon Length (AA):")
        exon_length_AA = st.sidebar.number_input("", min_value=None, max_value=None, value=5, step=None, format=None, key=None)

        if fasta1 !="" and fasta2 !="":
            st.markdown("#### Results")
            #st.write("\n")
            #st.markdown("##### Unfiltered Alignment:")
            #st.write("\n")
            maped_tuple = map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(fasta1, fasta2, match, mismatch, open_gap_penalty, gap_extension_penalty,exon_length_AA)
            Alignment_preview = visualise_unfiltered_alignment_better(maped_tuple[0])
            alignment_final = filter_and_visualize_final_alignment(maped_tuple[0])
            #st.text(Alignment_preview)
            st.write("\n")
            st.write("\n")
            st.markdown("##### Final Alignment:")
            st.write("\n")
            st.text(visualise_alignment_dynamically(maped_tuple[5],maped_tuple[6],maped_tuple[4]))
            st.write("\n")
           # st.text(maped_tuple[1])
           # st.text(maped_tuple[2])
           # st.text(maped_tuple[3])
           # st.text(maped_tuple[4])
           # st.text(maped_tuple[5])
           # st.text(maped_tuple[6])
            st.write("\n")
            st.markdown("##### Dataframe preview:")
            generated_table = write_results_to_tsv_file(maped_tuple,'/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv')
            st.write("\n")
            st.write(generated_table[1])
            st.write("\n")
            st.markdown("##### Download table:")
            st.write("\n")
            st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv', '','Mapped_Isoform_Positions.tsv'), unsafe_allow_html=True)
            st.write("--------------------------")


    elif choice == 'Download pre-computed data':
        st.header("Pre-computed mapped isoforms of GRCh37")
        st.write("--------------------------")
        st.markdown("#### Refseq (4GB):")
        st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv','', 'Refseq_Isoforms.tsv'), unsafe_allow_html=True)
        st.markdown("#### Ensembl (4GB):")
        st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv','', 'Ensembl_Isoforms.tsv'), unsafe_allow_html=True)
        st.markdown("#### All ID's (8GB):")
        st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv','', 'All_Isoforms.tsv'), unsafe_allow_html=True)

    elif choice == 'About & Source Code':
        st.header("Details")
        st.write("--------------------------")
        st.markdown("#### About this code:")
        st.write("\n")
        st.markdown("##### Needleman-Wunsch Algorithm:")
        st.write("It is also sometimes referred to as the optimal matching algorithm and the global alignment technique. The Needleman–Wunsch algorithm is still widely used for optimal global alignment, particularly when the quality of the global alignment is of the utmost importance.")
        st.markdown("##### check_for_wrong_exon_alignments function :")
        st.write(" This function helps to identify falsely aligned elements (distinct exons) when globally aligning isoforms of a gene. The Needleman Wunsch algorithm randomly alignes fractions of non-identical exons since the optimization of the algorithm is to maximize matches (which makes sense with homologues but not with isoforms). This function discards such fractions of the alignment by rejecting exons shorter than the defined minimal length (in AA).")
        st.markdown("#### Functions:")
        code = '''
        def transform_uploaded_data_type_accordingly(file):
            'uploaded files can be different types of files. A transformation is needed to interpret the data correctly
            Type of input: FASTA, FA and TXT
            Output type: depends on the case'
            
        '''
        st.code(code, language='python')



#Execution


#Default Needleman- Wunsch Parameters:
match=2
mismatch= -1.75
open_gap_penalty= -1
gap_extension_penalty= 0


#Execution
if __name__ == '__main__':
    main()


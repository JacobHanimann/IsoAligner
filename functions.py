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


#new function ideas

def select_canonical_sequence(isoforms):
 'function that can be dynamically applied to a set of isoform sequences'

#Back-end code (biological aspects)

def visualise_alignment_dynamically(reference_sequence_list,isoform_sequence_list,AA_match_evalutation_list,sequence1='sequence1: ', sequence2='sequence2: '):
    '''Function that returns an alignment of two sequences according to the AA_match_evalutation list generated from
     the check_for_wrong_exon_alignments() function in a visually pleasing fashion
     Input: 3 lists
     Output: 1 String (formatted with whitespace and newline character)'''
    correct_match_character = "|"
    wrong_match_character = "x"
    alignment_character_list = [" " if score =="gap" else correct_match_character if score=="correct" else wrong_match_character for score in AA_match_evalutation_list]
    whitespace = len(sequence1)*" "
    output_alignment_string= sequence1+''.join(reference_sequence_list)+'\n'+whitespace+''.join(alignment_character_list)+"\n"+sequence2+''.join(isoform_sequence_list)
    return output_alignment_string

def display_Needleman_Wunsch_Parameters_Sidebar():
    st.sidebar.markdown("### Function Parameters")
    st.sidebar.write("\n")
    st.sidebar.markdown("#### Minimal Exon Length (AA):")
    exon_length_AA = st.sidebar.number_input("", min_value=None, max_value=None, value=5, step=None,
                                             format=None, key=None)
    st.sidebar.write("\n")
    st.sidebar.markdown("#### Needleman-Wunsch Algorithm:")
    st.sidebar.write("\n")
    match = st.sidebar.number_input("match:", min_value=None, max_value=None, value=1, step=None, format=None,
                                    key=None)
    mismatch = st.sidebar.number_input("mismatch:", min_value=None, max_value=None, value=-2, step=None,
                                       format=None, key=None)
    open_gap_penalty = st.sidebar.number_input("open gap penalty:", min_value=None, max_value=None, value=-1.75,
                                               step=None, format=None, key=None)
    gap_extension_penalty = st.sidebar.number_input("gap extension penalty:", min_value=None, max_value=None,
                                                    value=0,
                                                    step=None, format=None, key=None)
    st.sidebar.write("\n")


def display_alignment_for_one_gene_from_database(reference_transcript,list_of_gene_objects,index_of_gene,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA,ID_type='ENST'):
    '''
    :param reference_transcript:
    :param list_of_gene_objects:
    :param index_of_gene:
    :return:
    '''
    for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
        if getattr(transcript, ID_type) == reference_transcript:
            reference_protein_sequence = getattr(transcript, "protein_sequence")
            break
    for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
        if getattr(transcript, ID_type) == reference_transcript:
            continue
        #st.text('Alignment')
        #st.text(reference_protein_sequence)
        #st.text(transcript.protein_sequence)
        mapped_tuple = map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(reference_protein_sequence,transcript.protein_sequence,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
        st.text(visualise_alignment_dynamically(mapped_tuple[5],mapped_tuple[6],mapped_tuple[4]))


def split_elements_from_user_input_string(string):
    '''
    Function that separates gene names/ID from each other
    :param string:
    :return: list of elements
    '''
    if "\n" in string:
        list_of_elements = list(string.split('\n'))
    elif ","  in string:
        list_of_elements = list(string.split(','))
    else:
        list_of_elements = [string]

    return list_of_elements


def identify_IDs_from_user_text_input(string):
    '''
    Function that identifies which ID's the user typed in with regex. Returns a dict of ID_types which can be used to search through the database more efficiently
    :param list of elements:
    :return: dict of ID_types
    '''
    dict_of_IDs = {}
    list_of_elements = split_elements_from_user_input_string(string)
    for element in list_of_elements:
        if re.search('ENSG\d+\.\d+',element):
            dict_of_IDs[element]='ENSG_version'
        elif re.search('ENSG\d{11}',element):
            dict_of_IDs[element]='ENSG'

    return dict_of_IDs

def search_through_database_with_known_ID_Type(list_of_gene_objects,dict_of_IDs):
    '''
    Function that searches trough database with gettatribute()
    :param database_list, list_of_IDs
    :return: dictionary of indexes of each element
    '''
    dict_element_indexes = {}
    for element,ID in dict_of_IDs.items():
        found = False
        parent_class = True
        if ID in ['ENSG_version','ENST','ENST_version','ENSP','ENSP_version']: #list must be completed
            parent_class = False
        for index,gene in enumerate(list_of_gene_objects):
            if found:
                break
            if parent_class:
                if getattr(gene,ID) == element:
                    dict_element_indexes[element] = index
                    break
            else:
                if type(gene.protein_sequence_isoform_collection) ==list:
                    for protein_sequence in gene.protein_sequence_isoform_collection:
                        if getattr(protein_sequence, ID) == element:
                            dict_element_indexes[element] = index
                            found = True
                            break
        else:
            dict_element_indexes[element] = 'not found'
    return dict_element_indexes

def generate_dictionary_of_canonical_IDs_from_user_input_IDs(dict_element_indexes,list_of_gene_objects):
    '''
    function that returns a dictionary in which each ID has the canonical ID as a value. If user specified ID, then key=value, if not, the canonical sequence must be extracted from the gene objecc
    :param ID_dictionary:
    :param list_of_gene_objects:
    :return: dictionary of canonical ID's used as a default reference
    '''


def fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,index_of_gene,ID="ENST"):
    '''
    function to get IDs of the gene object to be displayed in a selectbox to choose a reference
    :param  list_of_gene_objects, index_of_gene, optional:ID_type:
    :return: list
    '''
    list_of_transcripts=[getattr(sequence,ID) for sequence in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection]
    return list_of_transcripts


def transform_uploaded_data_type_accordingly(file):
    '''uploaded files can be different types of files. A transformation is needed to interpret the data correctly
    Type of input: FASTA, FA and TXT
    Output type: depends on the case'''

def extract_only_AA_of_Fasta_file(fasta_file):
    '''
    using Regex to match pattern of an AA sequence and creating a list in which each AA is an element
    Input type: String
    Output type: String
    '''
    without_newline = re.sub('\n', '', fasta_file)
    without_newline_and_whitespice = re.sub('\s', '', without_newline)
    sequence_of_AA_acronym = '[A-Z]+'
    minimal_length_of_AA_seq = '[A-Z]{10}'
    raw_AA_seq_list= re.findall(sequence_of_AA_acronym + minimal_length_of_AA_seq, without_newline_and_whitespice)
    if len(raw_AA_seq_list) >=1:
        return raw_AA_seq_list[0] #string
    else:
        return None


def map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(fasta1,fasta2,match,mismatch,gap_penalty,gap_extension_penalty,exon_length_AA):
    'maps FMI AA on COSMIC AA and creates list of AAposition and gaps'
    print('Mapping transcripts with Needleman-Wunsch...')
    clean_reference_fasta = extract_only_AA_of_Fasta_file(fasta1)
    alignments = pairwise2.align.globalms(clean_reference_fasta,extract_only_AA_of_Fasta_file(fasta2), match, mismatch, gap_penalty, gap_extension_penalty,one_alignment_only=True, penalize_extend_when_opening=True, penalize_end_gaps=False)
    alignment_COSMIC_fasta = list(alignments[0][0])
    alignment_isoform_fasta=list(alignments[0][1])
    isoform_pattern_check= check_for_wrong_exon_alignments(alignment_COSMIC_fasta,alignment_isoform_fasta,exon_length_AA)
    reference_position_list=[]
    isoform_positions_list=[]
    aminoacids=[]
    position_reference=1
    position_isoform=1
    for i in range(0,len(alignment_COSMIC_fasta)):
        if isoform_pattern_check[i]!='wrong':
            if alignment_COSMIC_fasta[i]==alignment_isoform_fasta[i]:
                aminoacids.append(alignment_COSMIC_fasta[i])
                reference_position_list.append(position_reference)
                isoform_positions_list.append(position_isoform)
                position_reference+=1
                position_isoform+=1
            if alignment_COSMIC_fasta[i]=='-' and alignment_isoform_fasta[i]=='-': #does it even happen?
                continue
            if alignment_COSMIC_fasta[i]!='-' and alignment_isoform_fasta[i]=='-':
                position_reference+=1
            if alignment_COSMIC_fasta[i]=='-' and alignment_isoform_fasta[i]!='-':
                position_isoform+=1
        else:
            position_reference += 1
            position_isoform += 1
    return(format_alignment(*alignments[0],full_sequences=True),aminoacids,reference_position_list,isoform_positions_list,isoform_pattern_check,alignment_COSMIC_fasta,alignment_isoform_fasta)


#This function also needs a sanity check
def check_for_wrong_exon_alignments(ref,isoform,exon_length_AA):     #Has to be adjusted so the minimal exon length can be varied with tool (more neighbours have to be counted)
    '''
    This function helps to identify falsely aligned elements (distinct exons) when globally aligning isoforms of a gene.
    The Needleman Wunsch algorithm also randomly alignes fractions of different exons but they do not represent the same aminoacid.
    since the optimization of the algorithm is to maximize matches (which makes sense with homologues but not with isoforms)
    :param ref: aligned sequence in form of a list
    :param isoform: aligned sequence in form of a list
    :return: list which categories each elements alignment in 'correct,wrong,gap'
    To do: Function could be written in fewer lines: few things are copied like category gap classification and appending the category..could be all shortened
    '''
    isoform_check=[]
    for index in range(0,len(ref)):
        score=0 #score which determines if the position is fullfills the minimal exon length parameter
        gap=False
        # start of array
        if index <= exon_length_AA:
            if ref[index] != isoform[index]:
                category = "gap"
                gap=True
            else: #same Aminoacid
                score +=1
                for sidestep in range(1,exon_length_AA): #check how many neighbours there are
                    if ref[index + sidestep] == isoform[index + sidestep]:
                        score +=1
                    else:
                        break
            if score >= exon_length_AA and gap!=True:
                category = 'correct'
            elif score <= exon_length_AA and gap!=True:
                category = 'wrong'
            isoform_check.append(category)
            continue
        #end of array
        if len(ref)-index <=exon_length_AA :
            if ref[index] != isoform[index]:
                category = "gap"
                gap = True
            else:  # same Aminoacid
                score += 1
                for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                    if ref[index - sidestep] == isoform[index - sidestep]:
                        score += 1
                    else:
                        break
            if score >= exon_length_AA and gap!=True:
                category = 'correct'
            elif score <= exon_length_AA and gap!=True:
                category = 'wrong'
            isoform_check.append(category)
            continue
        #middle of array, checks for matches both sided
        if ref[index]!=isoform[index]:
            category="gap"
            gap = True
        else:  # same Aminoacid
            score += 1
            for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                if ref[index + sidestep] == isoform[index + sidestep]:
                    score += 1
                else:
                    break
            for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                if ref[index - sidestep] == isoform[index - sidestep]:
                    score += 1
                else:
                    break
        if score >= exon_length_AA and gap!=True:
            category='correct'
        elif score <= exon_length_AA and gap!=True:
            category='wrong'
        isoform_check.append(category)
    return isoform_check


#This function also needs a sanity check
def check_for_wrong_exon_alignments_second(ref,isoform,exon_length_AA):     #Has to be adjusted so the minimal exon length can be varied with tool (more neighbours have to be counted)
    '''
    This function helps to identify falsely aligned elements (distinct exons) when globally aligning isoforms of a gene.
    The Needleman Wunsch algorithm also randomly alignes fractions of different exons but they do not represent the same aminoacid.
    since the optimization of the algorithm is to maximize matches (which makes sense with homologues but not with isoforms)
    :param ref: aligned sequence in form of a list
    :param isoform: aligned sequence in form of a list
    :return: list which categories each elements alignment in 'correct,wrong,gap'
    '''
    isoform_check=[]
    for index in range(0,len(ref)):
        score=0 #score which determines if the position is part of an exon which is at least 5 elements long
        gap=False
        #at the start and end of the array an element has just neighbours on one side
        if index <=4: #start of array
            if ref[index] != isoform[index]:
                category = "gap"
                gap=True
            else: score +=1
            if ref[index + 1] == isoform[index + 1]:
                score += 1
                if ref[index + 2] == isoform[index + 2]:
                    score += 1
                    if ref[index + 3] == isoform[index + 3]:
                        score += 1
                        if ref[index + 4] == isoform[index + 4]:
                            score += 1
                            if ref[index + 5] == isoform[index + 5]:
                                score += 1
                                if ref[index + 6] == isoform[index + 6]:
                                    score += 1
            if score >= exon_length_AA and gap!=True:
                category = 'correct'
            elif score <= exon_length_AA and gap!=True:
                category = 'wrong'
            isoform_check.append(category)
            continue

        #end of array
        if len(ref)-index <=4 :
            if ref[index] != isoform[index]:
                category = "gap"
                gap = True
            else:
                score += 1
            if ref[index - 1] == isoform[index - 1]:
                score += 1
                if ref[index - 2] == isoform[index - 2]:
                    score += 1
                    if ref[index - 3] == isoform[index - 3]:
                        score += 1
                        if ref[index - 4] == isoform[index - 4]:
                            score += 1
                            if ref[index - 5] == isoform[index - 5]:
                                score += 1
                                if ref[index - 6] == isoform[index - 6]:
                                    score += 1
            if score >= exon_length_AA and gap!=True:
                category = 'correct'
            elif score <= exon_length_AA and gap!=True:
                category = 'wrong'
            isoform_check.append(category)
            continue

        #middle of array, checks for matches both sided
        if ref[index]!=isoform[index]:
            category="gap"
            gap = True
        else:
            score += 1
        if ref[index + 1] == isoform[index + 1]: #neighbours to the right
            score += 1
            if ref[index + 2] == isoform[index + 2]:
                score += 1
                if ref[index + 3] == isoform[index + 3]:
                    score += 1
                    if ref[index + 4] == isoform[index + 4]:
                        score += 1
                        if ref[index + 5] == isoform[index + 5]:
                            score += 1
                            if ref[index + 6] == isoform[index + 6]:
                                score += 1
        if ref[index - 1] == isoform[index - 1]:  #neighbours to the left
            score += 1
            if ref[index - 2] == isoform[index - 2]:
                score += 1
                if ref[index - 3] == isoform[index - 3]:
                    score += 1
                    if ref[index - 4] == isoform[index - 4]:
                        score += 1
                        if ref[index - 5] == isoform[index - 5]:
                            score += 1
                            if ref[index - 6] == isoform[index - 6]:
                                score += 1
        if score >= exon_length_AA and gap!=True:
            category='correct'
        elif score <= exon_length_AA and gap!=True:
            category='wrong'
        isoform_check.append(category)
    return isoform_check


def write_results_to_tsv_file(mapped_tuple,file_location): #has to be redesign for different inputs
    'to be written'
    print('Writing results to csv file...')
    with open(file_location, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t')
        tsv_writer.writerow(['AA','ReferencePos','IsoformPos'])
        for indexiterator in range(0,len(mapped_tuple[1])):
            tsv_writer.writerow([mapped_tuple[1][indexiterator],mapped_tuple[2][indexiterator],mapped_tuple[3][indexiterator]])
    df = pd.read_csv(file_location, sep='\t')
    return (file_location,df)




#classes

class Gene:
    def __init__(self, ENSG, ensembl_gene_symbol,refseq_gene_ID=None, HGNC=None, HGNC_gene_symbol=None, previous_symbols=None, alias_symbols=None, protein_sequence_isoform_collection=None, canonical_default=None, average_exon_length=None, uniprot_ID=None):
        self.ENSG = ENSG
        self.ensembl_gene_symbol = ensembl_gene_symbol
        self.refseq_gene_ID = refseq_gene_ID
        self.HGNC = HGNC
        self.HGNC_gene_symbol = HGNC_gene_symbol
        self.previous_symbols = previous_symbols
        self.alias_symbols = alias_symbols
        self.protein_sequence_isoform_collection = protein_sequence_isoform_collection
        self.canonical_default = canonical_default
        self.average_exon_length= average_exon_length
        self.uniprot_ID = uniprot_ID


class protein_sequence:
    def __init__(self,gene_name, protein_sequence, ENSG=None, ENSG_version=None, ENST=None, ENST_version=None, ENSP=None,
                ENSP_version=None, refseq_rna=None, refseq_protein=None, uniprot_accession=None, uniprot_uniparc=None, uniprot_isoform=None):
        self.gene_name= gene_name #maybe unnecessary
        self.protein_sequence = protein_sequence
        self.ENSG = ENSG
        self.ENSG_version = ENSG_version
        self.ENST = ENST
        self.ENST_version = ENST_version
        self.ENSP = ENSP
        self.ENSP_version = ENSP_version
        self.refseq_rna = refseq_rna
        self.refseq_protein = refseq_protein
        self.uniprot_accession = uniprot_accession
        self.uniprot_uniparc = uniprot_uniparc
        self.uniprot_isoform = uniprot_isoform
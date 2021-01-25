import pandas as pd
import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np
from PIL import Image
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
from collections import Counter




#Functions for front-end streamlit components:

def get_binary_file_downloader_html(bin_file, file_label='File', name_button='download'):
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">{name_button} {file_label}</a>'
    return href


def get_table_download_link(df,name_of_file = "dataframe.tsv" ,sep='\t'):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=False, sep=sep)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}">'+name_of_file+'</a>'
    return href



#new function ideas

def select_canonical_sequence(isoforms):
 'function that can be dynamically applied to a set of isoform sequences'

#Back-end code (biological aspects)

def visualise_alignment_dynamically(reference_sequence_list,isoform_sequence_list,AA_match_evalutation_list,percentage_reference,percentage_isoform,sequence1='sequence1', sequence2='sequence2',):
    '''Function that returns an alignment of two sequences according to the AA_match_evalutation list generated from
     the check_for_wrong_exon_alignments() function in a visually pleasing fashion
     Input: 3 lists
     Output: 1 String (formatted with whitespace and newline character)'''
    correct_match_character = "|"
    wrong_match_character = "x"
    alignment_character_list = [" " if score =="gap" else correct_match_character if score=="correct" else wrong_match_character for score in AA_match_evalutation_list]
    whitespace = len(sequence1)*" "+"         "
    output_alignment_string= sequence1+":("+str(round(100*percentage_reference,1))+'%) '+''.join(reference_sequence_list)+'\n'+whitespace+''.join(alignment_character_list)+"\n"+sequence2+":("+str(round(100*percentage_isoform,1))+'%) '+''.join(isoform_sequence_list)
    return output_alignment_string


def display_alignment_for_one_gene_from_database(reference_transcript,list_of_gene_objects,index_of_gene,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA,ID_type='ENSP'):
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
    transcript_number = 1
    for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
        if getattr(transcript, ID_type) == reference_transcript:
            continue
        isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(reference_protein_sequence,transcript.protein_sequence,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)[4:7]
        percentage_reference, percentage_isoform = calculate_percentage_of_mapped_positions(isoform_pattern_check,reference_protein_sequence,transcript.protein_sequence)
        st.write('Alignment '+str(transcript_number))
        st.text(visualise_alignment_dynamically(alignment_reference_fasta,alignment_isoform_fasta,isoform_pattern_check,percentage_reference,percentage_isoform,sequence1 =reference_transcript,sequence2= transcript.ENST))
        st.text('\n')
        transcript_number +=1

def calculate_percentage_of_mapped_positions(isoform_check, reference_protein_sequence,isoform_protein_sequence):
    '''
    :return: percentage of mapped positions
    '''
    positions_total_reference = len(reference_protein_sequence)
    positions_total_isoform = len(isoform_protein_sequence)
    counter_dict = Counter(isoform_check)
    matches= counter_dict['correct']
    percentage_reference = matches / positions_total_reference
    percentage_isoform = matches / positions_total_isoform
    return percentage_reference,percentage_isoform



def split_elements_from_user_input_string(string):
    '''
    Function that separates gene names/ID from each other
    :param string:
    :return: list of elements
    '''
    string = re.sub(' ','',string)
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
    comment: should be updated for all kinds of ID in database
    '''
    dict_of_IDs = {}
    list_of_elements = split_elements_from_user_input_string(string)
    for element in list_of_elements:
        #ensembl
        if re.search('ENSG\d+\.\d+',element):
            dict_of_IDs[element]='ENSG_version'
            continue
        elif re.search('ENSG\d{11}',element):
            dict_of_IDs[element]='ENSG'
            continue
        if re.search('ENST\d+\.\d+',element):
            dict_of_IDs[element]='ENST_version'
            continue
        elif re.search('ENST\d{11}', element):
            dict_of_IDs[element] = 'ENST'
            continue
        if re.search('ENSP\d+\.\d+', element):
            dict_of_IDs[element] = 'ENSP_version'
            continue
        elif re.search('ENSP\d{11}', element):
            dict_of_IDs[element] = 'ENSP'
            continue
        #uniprot
        if re.search('[OPQ][0-9][0-9A-Z]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z,0-9]{2}[0-9]|[A-N,R-Z][0-9][A-Z][A-Z,0-9]{2}[0-9][A-Z][A-Z,0-9]{2}[0-9]', element):
            dict_of_IDs[element] = 'uniprot_accession'
            continue
        if re.search('UPI[0-9A-F]+',element):
            dict_of_IDs[element] = 'uniprot_uniparc'
            continue
        #refseq
        if re.search('NM_\d+\.\d+', element):
            dict_of_IDs[element] = 'refseq_rna_version'
            continue
        elif re.search('NM_\d+', element):
            dict_of_IDs[element] = 'refseq_rna'
            continue
        if re.search('NP_\d+\.\d+', element):
            dict_of_IDs[element] = 'refseq_prot_version'
            continue
        elif re.search('NP_\d+', element):
            dict_of_IDs[element] = 'refseq_prot'
            continue

        #if no ID's were found, the string is probably a gene name
        dict_of_IDs[element]= 'gene_name'

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
        if ID in ['ENSG_version','ENST','ENST_version','ENSP','ENSP_version', 'refseq_rna', 'refseq_protein','uniprot_accession','uniprot_uniparc','uniprot_isoform']: #list must be completed
            parent_class = False
        for index,gene in enumerate(list_of_gene_objects):
            if found:
                break
            if parent_class:
                if ID =="gene_name":
                    if getattr(gene, "ensembl_gene_symbol") == element:
                        dict_element_indexes[element] = index
                        break
                    if getattr(gene, "HGNC_gene_symbol") == element:
                        dict_element_indexes[element] = index
                        break
                    if type(gene.previous_symbols) == list: #line can be deleted since all these attributes should be lists
                        if element in getattr(gene, "previous_symbols"):
                            dict_element_indexes[element] = index
                            break
                    if type(gene.alias_symbols)==list:
                        if element in getattr(gene, "alias_symbols"):
                            dict_element_indexes[element] = index
                            break
                else:
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


def fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,index_of_gene,ID="ENSP"):
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
    alignment_reference_fasta = list(alignments[0][0])
    alignment_isoform_fasta=list(alignments[0][1])
    isoform_pattern_check= check_for_wrong_exon_alignments(alignment_reference_fasta, alignment_isoform_fasta, exon_length_AA)
    reference_position_list=[]
    isoform_positions_list=[]
    aminoacids=[]
    position_reference=1
    position_isoform=1
    for i in range(0, len(alignment_reference_fasta)):
        if isoform_pattern_check[i]!='wrong':
            if alignment_reference_fasta[i]==alignment_isoform_fasta[i]:
                aminoacids.append(alignment_reference_fasta[i])
                reference_position_list.append(position_reference)
                isoform_positions_list.append(position_isoform)
                position_reference+=1
                position_isoform+=1
            if alignment_reference_fasta[i]== '-' and alignment_isoform_fasta[i]== '-': #does it even happen?
                continue
            if alignment_reference_fasta[i]!= '-' and alignment_isoform_fasta[i]== '-':
                position_reference+=1
            if alignment_reference_fasta[i]== '-' and alignment_isoform_fasta[i]!= '-':
                position_isoform+=1
        else:
            position_reference += 1
            position_isoform += 1
    return(format_alignment(*alignments[0],full_sequences=True), aminoacids, reference_position_list, isoform_positions_list, isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta)


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


def sidebar_pop_up_parameters():
    '''
    sidebar pop up for the parameters of all functions which can be dynamically adjusted.
    '''
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

    return match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA



#def write_results_to_tsv_file(mapped_tuple,file_location): #has to be redesign for different inputs
#    'to be written'
#    print('Writing results to csv file...')
#    with open(file_location, 'w') as output_file:
#        tsv_writer = csv.writer(output_file, delimiter='\t')
#        tsv_writer.writerow(['AA','ReferencePos','IsoformPos'])
#        for indexiterator in range(0,len(mapped_tuple[1])):
#            tsv_writer.writerow([mapped_tuple[1][indexiterator],mapped_tuple[2][indexiterator],mapped_tuple[3][indexiterator]])
#    df = pd.read_csv(file_location, sep='\t')
#    return (file_location,df)



def create_pandas_dataframe_raw_aa_sequence(needleman_mapped):
    nested_list = [[needleman_mapped[1][indexiterator],needleman_mapped[2][indexiterator],needleman_mapped[3][indexiterator]] for indexiterator in range(0,len(needleman_mapped[1]))]
    df = pd.DataFrame(nested_list, columns=(['AA', 'sequence1', 'sequence2']))
    return df


def create_table_for_one_gene_object(chosen_reference,list_of_gene_objects, index_of_gene,chosen_columns, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA, ID_type="ENSP"):
    '''
    function to create a pandas dataframe for only one gene object
    :param chosen_reference:
    :param list_of_gene_objects:
    :param index_of_gene:
    :param chosen_columns:
    :return: pandas dataframe
    '''
    list_of_all_alignments = []
    for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
        if getattr(transcript, ID_type) == chosen_reference:
            reference_protein_sequence = getattr(transcript, "protein_sequence")
            break
    for index,transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
        if getattr(transcript, ID_type) == chosen_reference:
            index_reference_transcript = index
            continue
        aminoacids, reference_position_list, isoform_positions_list = map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)[1:4]

        def get_selected_columns_attributes_and_column_names(chosen_columns):
            '''
            function to get the column values selected for each row
            :param chosen_columns: from multiselectbox
            '''
            column_values = []
            column_names= []

            #most column names missing

            if "Gene name" in chosen_columns:
                column_values.append(list_of_gene_objects[index_of_gene].ensembl_gene_symbol)
                column_names.append("Gene_name")

            if "Ensembl Gene ID" in chosen_columns:
                column_values.append(transcript.ENSG)
                column_names.append("ENSG")

            if "Ensembl Transcript ID" in chosen_columns:
                column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[index_reference_transcript].ENST)
                column_names.append("Ref_transcript_ID")
                column_values.append(transcript.ENST)
                column_names.append("Iso_transcript_ID")

            if "Ensembl Protein ID" in chosen_columns:
                column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                         index_reference_transcript].ENSP)
                column_names.append("Ref_protein_ID")
                column_values.append(transcript.ENSP)
                column_names.append("Iso_protein_ID")

            column_names = column_names + ['AA','ReferencePos', 'IsoformPos']
            return column_values, column_names

        for indexiterator in range(0, len(aminoacids)):
            column_values,column_names = get_selected_columns_attributes_and_column_names(chosen_columns)
            positions = [aminoacids[indexiterator], reference_position_list[indexiterator], isoform_positions_list[indexiterator]]
            nested_list_alignment=column_values+positions
            list_of_all_alignments.append(nested_list_alignment)

    df = pd.DataFrame(list_of_all_alignments, columns=(column_names))

    return df


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
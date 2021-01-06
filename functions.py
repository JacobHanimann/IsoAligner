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

def visualise_alignment_dynamically(reference_sequence_list,isoform_sequence_list,AA_match_evalutation_list):
    '''Function that returns an alignment of two sequences according to the AA_match_evalutation list generated from
     the check_for_wrong_exon_alignments() function in a visually pleasing fashion
     Input: 3 lists
     Output: 1 String (formatted with whitespace and newline character)'''
    correct_match_character = "|"
    wrong_match_character = "x"
    alignment_character_list = [" " if score =="gap" else correct_match_character if score=="correct" else wrong_match_character for score in AA_match_evalutation_list]
    output_alignment_string= ''.join(reference_sequence_list)+'\n'+''.join(alignment_character_list)+"\n"+''.join(isoform_sequence_list)
    return output_alignment_string


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
    minimal_length_of_AA_seq = '[A-Z]{15}'
    raw_AA_seq= re.findall(sequence_of_AA_acronym + minimal_length_of_AA_seq, without_newline_and_whitespice)[0]
    return raw_AA_seq #string

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
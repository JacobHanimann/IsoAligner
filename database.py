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
from functions import * #import all functions from the functions.py file


#keep in mind that code also has to be descriptive to generate pre-computed offline data and not only for the dynamic stuff


class gene:
    def __init__(self, HGNC, gene_symbol, previous_symbols=None, alias_symbols=None, protein_sequence_isoform_collection=None, canonical_default=None):
        self.HGNC = HGNC
        self.gene_symbol = gene_symbol
        self.previous_symbols = previous_symbols
        self.alias_symbols = alias_symbols
        self.protein_sequence_isoform_collection = protein_sequence_isoform_collection
        self.canonical_default = canonical_default


def create_list_of_gene_objects(file_of_gene_names):
    '''make a list of gene objects
    input: textfile
    output: list of gene objects
    '''
    df = pd.read_csv(file_of_gene_names, sep='\t')
    list_of_gene_objects = [gene(df.loc[index,'HGNC'], gene_symbol = df.loc[index, 'approved_symbol'],previous_symbols = df.loc[index, 'previous_symbols'], alias_symbols = df.loc[index, 'alias_symbols'])for index in range(0,len(df))]
    for gene_object in list_of_gene_objects: #convert comma separated strings into elements of a list to facilitate a infrastructure which can be better searched through later (no need for regex later)
       if type(gene_object.previous_symbols) != float: #None values are type float
           if "," in  gene_object.previous_symbols:
            gene_object.previous_symbols = gene_object.previous_symbols.split(', ')
       if type(gene_object.alias_symbols) != float:
           if "," in gene_object.alias_symbols:
            gene_object.alias_symbols = gene_object.alias_symbols.split(', ')
    for gen in list_of_gene_objects:
        print(gen.HGNC, gen.gene_symbol,'first',gen.alias_symbols,'other', gen.previous_symbols)
    return list_of_gene_objects


class protein_sequence:
    def __init__(self,gene_name, protein_sequence, ENSG=None, ENSG_version=None, ENST=None, ENST_version=None, ENSP=None,
                ENSP_version=None, refseq_rna=None, refseq_protein=None, uniprot_accession=None, uniprot_uniparc=None, average_exon_length=None):
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
        self.average_exon_length= average_exon_length

def get_ensembl_fasta_sequences_and_IDs(file):
    '''What we need: reading in chunk and extracting AA sequence and ID's and then create a protein_sequence object
    try to be generic with regular expressions for fetching the ID's
    ouput: objects in a dictionary or a list'''


def get_refseq_fasta_sequences_and_IDs(file, list_of_objects):
    'also get refseq data'

    def check_if_IDs_can_be_mapped():
        'before creating new protein_sequence, check if it already exists'


def get_bio_IDs_with_regex(ID_type,string):
    'write generic functions to extract certain ID types from different databases'

    #Ensembl
    if ID_type=='ensembl_ensg':
        pattern = 'ENSG\d+{17}'
    elif ID_type == 'ensembl_ensg_version':
        pattern = 'ENSG\d+\.\d'
    elif ID_type == 'ensembl_enst':
         pattern = 'ENST\d+{17}'
    elif ID_type == 'ensembl_enst_version':
         pattern = 'ENST\d+\.\d'
    elif ID_type == 'ensembl_ensp':
         pattern = 'ENSP\d+{17}'
    elif ID_type == 'ensembl_ensp_version':
        pattern = 'ENSP\d+\.\d'

    #Refseq
    elif ID_type=='refseq_rna':
         pattern = 'NM_\d+'
    elif ID_type=='refseq_rna_version':
         pattern = 'NM_\d+\.\d+'
    elif ID_type=='refseq_prot':
         pattern = 'NP_\d+'
    elif ID_type=='refseq_prot_version':
        pattern = 'NP_\d+\.\d+'

    #Uniprot IDs
    elif ID_type == 'uniprot_accession':
         pattern = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
    elif ID_type == 'uniprot_uniparc':
        pattern = 'UPI[0-9A-F]+'

    match_list = re.findall(pattern,string)
    if not match_list: #if list is empty
        return 'not found'
    elif len(match_list)==1:
        return match_list[0]
    else:
        return match_list


def select_canonical_sequence(isoforms):
 'function that can be dynamically applied to a set of isoform sequences'

def save_all_data_in_pickle_style():
    'see if this method is fast enough with streamlit'


def map_isoforms_with_canonical_sequence():
    'I dont know if this makes sense'


def save_results_to_tsv_file(dictionary):
    'to be pre-computed values'


#Execution

list_of_gene_objects = create_list_of_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/HGNC_protein_coding.txt')

print(get_bio_IDs_with_regex('ensembl_enst', 'ENSG00000001036|ENSG00000001036.14|ENST00000002165|ENST00000002165.11|ENSP00000002165|ENSP00000002165.5|Q9BTY2|UPI0000073C10|FUCA2'))

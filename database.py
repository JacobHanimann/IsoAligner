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

#create a dictionary of genes (gene names in list and then make a search with if gene_name matches any element in the gene name list...
#all protein_Sequence objects in value of dictionary, then, use grouping to choose a canonical sequence or let the user select


def create_gene_dictionary(list_of_gene_names):
    'make a dictionary with all the alternative gene name per gene'

class gene:
    def __init__(self, HGNC, gene_symbol, previous_symbol, alias_symbol, protein_sequence_isoform_collection):
        self.HGNC = HGNC
        self.gene_symbol = gene_symbol
        self.previous_symbol = previous_symbol
        self.alias_symbol = alias_symbol
        self.protein_sequence_isoform_collection = protein_sequence_isoform_collection


class protein_sequence:
    def __init__(self,gene_name, protein_sequence, ENSG, ENSG_version, ENST, ENST_version, ENSP, ENSP_version, refseq_rna, refseq_protein, uniprot_accession, uniprot_uniparc, average_exon_length):
        self.gene_name= gene_name
        self.protein_sequence
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


def get_bio_IDs_with_regex(ID_type):
    'write generic functions to extract certain ID types from different databases'



def select_canonical_sequence(isoforms):
 'function that can be dynamically applied to a set of isoform sequences'

def select_and_map_isoforms_with_canonical_sequence():
    'I dont know if this makes sense'


def save_all_data_in_pickle_style():
    'see if this method is fast enough with streamlit'


def save_results_to_tsv_file(dictionary):
    'to be pre-computed values'
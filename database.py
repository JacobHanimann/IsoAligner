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


#two instances of searching for a sequence, first by gene name and then by ID

class protein_sequence:
    def __init__(self,gene_name, ENSG, ENSG_version, ENST, ENST_version, ENSP, ENSP_version, refseq_rna, refseq_protein, uniprot_accession, uniprot_uniparc, average_exon_length):
        self.gene_name= gene_name
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

def select_canonical_sequence(isoforms):
 'function that can be dynamically applied to a set of isoform sequences'


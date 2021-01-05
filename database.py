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

#Allgemeins Problem: an welem kriterium sell ich es gene definiere: an name? anere ID?

#create a dictionary of genes (gene names in list and then make a search with if gene_name matches any element in the gene name list...

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


#choose the canonical sequence

#generate the 3 three mapped lists and store them somewhere (only used when a gene name is entered or the user takes the canonical one

#what if the user wants a different reference protein sequence?
# I have to then compute it realtime probably
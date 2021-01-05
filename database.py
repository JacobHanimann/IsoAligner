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

#create a dictionary of references (ensembl, refseq, uniprot)     Problem: Es git überschnidige

#Allgemeins Problem: an welem kriterium sell ich es gene definiere: an name? anere ID?

#create a dictionary of genes

class protein_sequence:
    def __init__(self,gene_name, ENSG, ENSG_version, ENST, ENST_version, ENSP, ENSP_version, uniprot_accession, uniprot_uniparc):
        self.gene_name= gene_name
        self.ENSG = ENSG
        self.ENSG_version = ENSG_version

#choose the canonical sequence

#generate the 3 three mapped lists and store them somewhere (only used when a gene name is entered or the user takes the canonical one

#what if the user wants a different reference protein sequence?
# I have to then compute it realtime probably


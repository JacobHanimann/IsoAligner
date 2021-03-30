import pandas as pd
import re
from Gene import *
from Protein_isoform import *
from Alignment import *
from Extractions_BioIDs import *

class Exon_Information():
    pass

    @staticmethod
    def read_Ensembl_GRCh38_gtf_file(file):
        '''extract fasta files one by one and add them to the gene objects'''
        with open(file, "r") as f:
            expenses_txt = f.readlines()
        for line in expenses_txt:
            print(line)
        whole_txt = "".join(expenses_txt)
        splittext = re.split(">", whole_txt)
        fasta_count = 0
        matches = 0
        list_of_gene_objects = []
        for fasta in splittext[1:]:
            fasta_count += 1


#Execution

Exon_Information.read_Ensembl_GRCh38_gtf_file('/Users/jacob/Desktop/Isoform Mapper Webtool/Homo_sapiens.GRCh38_protein_coding.gtf')


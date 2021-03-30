import pandas as pd
import re
from Gene import *
from Protein_isoform import *
from Alignment import *
from Extractions_BioIDs import *
from Exon import *

class Exon_Information():
    pass

    @staticmethod
    def read_Ensembl_GRCh38_gtf_file_generate_nested_dict(file):
        '''extract fasta files one by one and add them to the gene objects'''
        with open(file, "r") as f:
            all_lines = f.readlines()
        gene_dict={}
        exon_dict={}
        list_of_nested_dict =[]
        for line in all_lines:
            Exon()



    @staticmethod
    def pick_exon_length_median_from_nested_dict():
        pass



def add_exon_median_to_gene_objects():
    pass

def add_exon_objects_to_protein_objects():
    pass


#Execution

Exon_Information.read_Ensembl_GRCh38_gtf_file('/Users/jacob/Desktop/Isoform Mapper Webtool/Homo_sapiens.GRCh38_protein_coding.gtf')

#nested dict as an output

#newdict with median of each ENSG

#assign new dict to gene object

#assign exon objects to protein isoforms if wanted
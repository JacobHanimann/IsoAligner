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
        #read file
        with open(file, "r") as f:
            all_lines = f.readlines()
        #organisation
        gene_dict={}
        print(gene_dict)
        for line in all_lines:
            splitted = re.split('\t',line)
            #extracting information
            ENSG = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_ensg',splitted[8])
            start_end = [int(splitted[3]),int(splitted[4])]
            exon_length = (start_end[1]-start_end[0]+1)/3 #Lösung finden
            print('exon length:', exon_length)
            ENSE = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_ense',splitted[8])
            if ENSE == 'not found':
                ENSE = None
            ENST = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_enst',splitted[8])
            exon_number_extract = re.findall('exon_number\s"\d"',splitted[8])
            if exon_number_extract:
               exon_number = re.findall('\d',exon_number_extract[0])[0]
            else:
                exon_number = None
            #creating Exon_objects
            exon_object = Exon(exon_start_end=start_end,exon_length=exon_length,ENSE=ENSE,ENST=ENST,exon_number=exon_number)
            #saving exon object in gene dict
            if ENSG in gene_dict:
                gene_dict[ENSG].append(exon_object)
            else:
                gene_dict[ENSG] = [exon_object]
        return gene_dict


    @staticmethod
    def pick_exon_length_median_from_nested_dict():
        pass


def add_exon_median_to_gene_objects():
    pass

def add_exon_objects_to_protein_objects():
    pass


#Execution

Exon_Information.read_Ensembl_GRCh38_gtf_file_generate_nested_dict('/Users/jacob/Desktop/Isoform Mapper Webtool/Homo_sapiens.GRCh38_protein_coding.gtf')

#nested dict as an output

#newdict with median of each ENSG

#assign new dict to gene object

#assign exon objects to protein isoforms if wanted
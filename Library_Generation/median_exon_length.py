import pandas as pd
import re
from Gene import *
from Protein_isoform import *
from Alignment import *
from Extractions_BioIDs import *
from Exon import *
import statistics
import pickle

class Exon_Information():
    pass

    @staticmethod
    def read_Ensembl_GRCh38_gtf_file_generate_nested_dict(file):
        '''extraction of exon infos of protein coding genes, generates exon objects and stores them in dictionary with ENSG IDs as keys'''
        #read file
        with open(file, "r") as f:
            all_lines = f.readlines()
        #organisation
        gene_dict={}
        exon_length = None
        ENSE= None
        cds_exon_number = None
        exon_number = None
        print('lines in total:', len(all_lines))
        for index, line in enumerate(all_lines):
            if index % 100000==0:
                print(100*round(index/len(all_lines),2),'%')
            splitted = re.split('\t',line)

            if splitted[2]=='exon':
            #extracting information
                ENSG = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_ensg',splitted[8])
                start_end = [int(splitted[3]),int(splitted[4])]
                ENSE = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_ense',splitted[8])
                if ENSE == 'not found':
                    ENSE = None
                ENST = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_enst',splitted[8])
                exon_number_extract = re.findall('exon_number\s"\d"',splitted[8])
                if exon_number_extract:
                    exon_number = re.findall('\d',exon_number_extract[0])[0]
                else:
                    exon_number = None

            if splitted[2]=="CDS":
                cds_start_end = [int(splitted[3]), int(splitted[4])]
                frame = int(splitted[7])
                #print('new cds')
                #print(frame)
                #print(cds_start_end)
                exon_length = (cds_start_end[1]-cds_start_end[0]+1-frame)/3 #stimmt immer noch nicht, mit Abdullah besprechen
                #print('exon length:', exon_length)
                exon_number_extract = re.findall('exon_number\s"\d"', splitted[8])
                if exon_number_extract:
                    cds_exon_number = re.findall('\d', exon_number_extract[0])[0]
                else:
                    cds_exon_number = None


            #creating Exon_objects
            if exon_length!=None and ENSE!=None and cds_exon_number!=None and cds_exon_number==exon_number:
                exon_object = Exon(exon_start_end=start_end,exon_length=exon_length,ENSE=ENSE,ENST=ENST,exon_number=exon_number)
                #saving exon object in gene dict
                if ENSG in gene_dict:
                    gene_dict[ENSG].append(exon_object)
                else:
                    gene_dict[ENSG] = [exon_object]
                exon_length=None
                ENSE=None
                cds_exon_number = None
                exon_number = None
                #print(splitted)
                #print(exon_object.__dict__)
        return gene_dict


    @staticmethod
    def pick_exon_length_median_from_nested_dict(gene_dict):
        '''
        functions that takes the median of the exon length per gene and stores it in a dict
        :param gene_dict: dictionary generated from read ensembl gtf file
        :return: dict
        '''
        gene_dict_only_median = {}
        for gene,exon_list in gene_dict.items():
            exon_lengths = [exon.exon_length_in_AA for exon in exon_list]
            median_of_gene = statistics.median(exon_lengths)
            gene_dict_only_median[gene] = median_of_gene

        return gene_dict_only_median


    @staticmethod
    def add_exon_median_to_gene_objects(list_of_gene_objects,gene_dict_median):
        print('total length:',len(gene_dict_median))
        index = 0
        for ENSG, exon_length in gene_dict_median.items():
            index += 1
            if index % 1000 == 0:
                print(100 * round(index / len(gene_dict_median), 2), '%')
            for gene in list_of_gene_objects:
                if ENSG==gene.ENSG:
                    gene.median_exon_length = round(exon_length) #has to be adjusted
                    break


    @staticmethod
    def add_exon_objects_to_protein_objects(list_of_gene_objects,gene_dict):
        index = 0
        for ENSG, exon_list in gene_dict.items():
            index += 1
            if index % 1000 == 0:
                print(100 * round(index / len(gene_dict), 2), '%')
            for gene in list_of_gene_objects:
                if ENSG == gene.ENSG:
                    for exon in exon_list:
                        for isoform in gene.protein_sequence_isoform_collection:
                            if exon.ENST==isoform.ENST:
                                if isoform.collection_of_exons ==None:
                                    isoform.collection_of_exons = [exon]
                                else:
                                    isoform.collection_of_exons.append(exon)
                                break
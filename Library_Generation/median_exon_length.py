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
        print('lines in total:', len(all_lines))
        for index, line in enumerate(all_lines):
            if index % 100000==0:
                print(100*round(index/len(all_lines),2),'%')
            splitted = re.split('\t',line)
            #print(splitted)
            #extracting information
            ENSG = Get_Bio_ID.get_bio_IDs_with_regex('ensembl_ensg',splitted[8])
            start_end = [int(splitted[3]),int(splitted[4])]
            exon_length = (start_end[1]-start_end[0]+1)/3 #Lösung finden
            #print('exon length:', exon_length)
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
                    gene.median_exon_length = exon_length
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
                                print('zugewiesen')
                                break


#Execution
gene_dict = Exon_Information.read_Ensembl_GRCh38_gtf_file_generate_nested_dict('/Users/jacob/Desktop/Isoform Mapper Webtool/Homo_sapiens_GRCh38_exons.gtf')

#print('Pickling genes dict')
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/genes_dict", "wb") as fp:  # Pickling
#    pickle.dump(gene_dict, fp)

genes_dict_median = Exon_Information.pick_exon_length_median_from_nested_dict(gene_dict)


with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_25_march_fourth.txt", "rb") as fp:  # Pickling
    list_of_gene_objects = pickle.load(fp)

#Exon_Information.add_exon_median_to_gene_objects(list_of_gene_objects,genes_dict_median)
Exon_Information.add_exon_objects_to_protein_objects(list_of_gene_objects,gene_dict)

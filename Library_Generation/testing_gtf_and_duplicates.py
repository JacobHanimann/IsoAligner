from Ensembl import *
from HGNC import *
from Biomart_tables import *
from Refseq import *
from Uniprot import *
from Validation_of_library import *
from minimal_exon_length import *

date = '27th_may'




print('Pickling list of gene objects and saving file...')
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_fifth.txt", "rb") as fp:  # Pickling
    list_of_gene_objects = pickle.load(fp)


def check_if_there_are_AA_seq_duplicates(list_of_gene_objects):
    '''
    check out if there were IDs and Seq that are the same but escaped the match
    :param list_of_gene_objects:
    :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
    '''

    def duplicates(lst, item):
        return [i for i, x in enumerate(lst) if x == item]

    genes_without_AA_seq = 0
    duplicates_number = 0
    genes_without_duplicates = 0
    genes_with_more_than_one_duplicate = 0
    redundant_sequences = 0
    duplicate_genes_dict = dict()
    list_of_all_seq= []
    for index, gene in enumerate(list_of_gene_objects):
        if type(gene.protein_sequence_isoform_collection) == list:
            List = [sequence.protein_sequence for sequence in gene.protein_sequence_isoform_collection]
            list_of_all_seq = list_of_all_seq + List
        print(index)
    duplicates_dict = dict((x, duplicates(list_of_all_seq, x)) for x in set(list_of_all_seq) if list_of_all_seq.count(x) > 1)


    print(duplicates_dict)
    print(len(duplicates_dict))



check_if_there_are_AA_seq_duplicates(list_of_gene_objects)
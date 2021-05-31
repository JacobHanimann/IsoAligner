from Ensembl import *
from HGNC import *
from Biomart_tables import *
from Refseq import *
from Uniprot import *
from minimal_exon_length import *

date = '27th_may'




print('Pickling list of gene objects and saving file...')
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_fifth.txt", "rb") as fp:  # Pickling
    list_of_gene_objects = pickle.load(fp)


def check_if_there_are_AA_seq_duplicates_over_all_genes(list_of_gene_objects):
    '''
    check out if there were IDs and Seq that are the same but escaped the match
    :param list_of_gene_objects:
    :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
    '''

    def duplicates(lst, item):
        return [i for i, x in enumerate(lst) if x == item]
    list_of_all_names= []
    for index, gene in enumerate(list_of_gene_objects):
        list_of_all_names.append(gene.uniprot_name_ID)
    duplicates_dict = dict((x, duplicates(list_of_all_names, x)) for x in set(list_of_all_names) if list_of_all_names.count(x) > 1)

    print(duplicates_dict)
    print('gene object with same name:',len(duplicates_dict))
    return duplicates_dict


check_if_there_are_AA_seq_duplicates_over_all_genes(list_of_gene_objects)
import streamlit as st

class Statistics:
    pass

    @staticmethod
    def list_of_gene_objects_statistics(list_of_gene_objects):
        '''
        function that returns statistics about the whole library
        :param list_of_gene_objects:
        :return: dictionary
        '''
        total_number_of_genes = len(list_of_gene_objects)
        total_number_of_isoforms = 0
        genes_without_isoforms = 0
        for gene in list_of_gene_objects:
            if type(gene.protein_sequence_isoform_collection) == list:
                total_number_of_isoforms = total_number_of_isoforms + len(gene.protein_sequence_isoform_collection)
                ensembl_ids = 0  # function that counts ensembl ids per gene object as a method of the gene class
                refseq_ids = 0
                uniprot_ids = 0
            else:
                genes_without_isoforms += 1
        return total_number_of_genes, total_number_of_isoforms, genes_without_isoforms


    @staticmethod
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
        for index, gene in enumerate(list_of_gene_objects):
            if type(gene.protein_sequence_isoform_collection) == list:
                List = [sequence.protein_sequence for sequence in gene.protein_sequence_isoform_collection]
                duplicates_dict = dict((x, duplicates(List, x)) for x in set(List) if List.count(x) > 1)
                if len(duplicates_dict) != 0:
                    if list(duplicates_dict.keys())[0] != None:
                        duplicate_genes_dict[index] = duplicates_dict
                        duplicates_number += 1
                        print('new gene now: ', gene.ensembl_gene_symbol)
                        print(duplicates_dict)
                        for sequence, objects in duplicates_dict.items():
                            redundant_sequences = redundant_sequences + len(objects)
                    if len(duplicates_dict) > 1:
                        genes_with_more_than_one_duplicate += 1
                else:
                    genes_without_duplicates += 1
            else:
                genes_without_AA_seq += 1

        return duplicates_number, genes_without_duplicates, redundant_sequences, genes_with_more_than_one_duplicate


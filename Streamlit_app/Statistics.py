import streamlit as st
from Protein_isoform import *
import statistics
from Gene import *

class Statistics:
    pass

    @staticmethod
    def list_of_gene_objects_statistics(list_of_gene_objects):
        '''
        function that returns statistics about the whole library
        :param list_of_gene_objects:
        :return: dictionary
        '''
        number_of_genes = 0
        more_than_two_isoforms = 0
        total_number_of_isoforms = 0
        number_of_isoforms_two = 0
        genes_without_mme = 0
        IDs_in_total = 0
        number_of_IDs_two = 0
        minimal_exon_lengths = []
        isoform_attributes = Protein_isoform.list_of_attributes()
        gene_attributes = Gene.list_of_attributes()
        gene_attributes.remove('protein_sequence_isoform_collection')
        for gene in list_of_gene_objects:
            if gene.minimal_exon_length != None:
                number_of_genes +=1
                if len(gene.protein_sequence_isoform_collection) >= 2:
                     more_than_two_isoforms += 1
                     number_of_isoforms_two = number_of_isoforms_two + len(gene.protein_sequence_isoform_collection)
                     for attribute in gene_attributes:
                        if getattr(gene, attribute) != None:
                            number_of_IDs_two += 1
                     for sequence in gene.protein_sequence_isoform_collection:
                        for attribute in isoform_attributes:
                            if getattr(sequence,attribute)!=None:
                                number_of_IDs_two +=1
                minimal_exon_lengths.append(gene.minimal_exon_length)
                for attribute in gene_attributes:
                    if getattr(gene, attribute) != None:
                        IDs_in_total += 1
                if type(gene.protein_sequence_isoform_collection) == list:
                    total_number_of_isoforms = total_number_of_isoforms + len(gene.protein_sequence_isoform_collection)
                    for sequence in gene.protein_sequence_isoform_collection:
                        for attribute in isoform_attributes:
                            if getattr(sequence,attribute)!=None:
                                IDs_in_total +=1
            else:
                genes_without_mme += 1
        minimal_exon_lengths = [exon for exon in minimal_exon_lengths if exon >=0]
        median_exon = statistics.median(minimal_exon_lengths)
        return number_of_genes, more_than_two_isoforms, total_number_of_isoforms, genes_without_mme, IDs_in_total, len(minimal_exon_lengths), round(median_exon), number_of_isoforms_two,number_of_IDs_two


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
                        for sequence, objects in duplicates_dict.items():
                            redundant_sequences = redundant_sequences + len(objects)
                    if len(duplicates_dict) > 1:
                        genes_with_more_than_one_duplicate += 1
                else:
                    genes_without_duplicates += 1
            else:
                genes_without_AA_seq += 1

        return duplicates_number, genes_without_duplicates, redundant_sequences, genes_with_more_than_one_duplicate

    @staticmethod
    def isoform_form_check_stats(isoform_check_list):
        correct_aa = 0
        false_aa = 0
        mismatch_aa =0
        for type in isoform_check_list:
            if type == 'gap':
                continue
            if type == 'correct':
                correct_aa += 1
            if type == 'wrong':
                false_aa += 1
            if type =='mismatch':
                mismatch_aa +=1
        return correct_aa, false_aa,mismatch_aa


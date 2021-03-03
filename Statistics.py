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
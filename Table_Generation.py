import pandas as pd
from Alignment import *
import streamlit as st

class Table_Generation:
    pass

    @staticmethod
    def create_pandas_dataframe_raw_aa_sequence(needleman_mapped):
        '''input: 3 lists generated from the map_Needleman Wunsch function
        output: pandas dataframe'''
        nested_list = [
            [needleman_mapped[1][indexiterator], needleman_mapped[2][indexiterator], needleman_mapped[3][indexiterator]]
            for indexiterator in range(0, len(needleman_mapped[1]))]
        df = pd.DataFrame(nested_list, columns=(['AA', 'sequence1', 'sequence2']))
        return df


    @staticmethod
    def create_table_for_one_gene_object(index_reference_transcript, list_of_gene_objects, index_of_gene, chosen_columns, match,
                                         mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA,
                                         ID_type="ENSP", one_ID=True):
        '''
        one_ID: function returns pandas dataframe directly of alignments
        multiple IDs: function returns with list of the computed alignments which is then further used in create_table_for_dict_of_gene_objects function

        :param chosen_reference: index
        :param list_of_gene_objects: (library)
        :param index_of_gene:
        :param chosen_columns: (choosed from selectbox)
        :return: pandas dataframe or lists of computed alignment
        '''
        list_of_all_alignments = []

        #select reference protein sequence
        reference_protein_sequence = list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[index_reference_transcript].protein_sequence  # chosen_reference is an index (table generation for multiple ID's)


        #create alignment for each alternative isoform
        for index, transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
            if index == index_reference_transcript: #do not align the reference transcript with itself
                continue
            aminoacids, reference_position_list, isoform_positions_list = Alignment.map_AA_Needleman_Wunsch_with_exon_check(
                reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)[1:4]


            def get_selected_columns_attributes_and_column_names(chosen_columns):
                '''
                function to get the column values selected for each row
                :param chosen_columns: from multiselectbox
                '''
                column_values = []
                column_names = []

                # most column names missing

                if "Gene name" in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].ensembl_gene_symbol)
                    column_names.append("Gene_name")

                if "Ensembl Gene ID (ENSG)" in chosen_columns:
                    column_values.append(transcript.ENSG)
                    column_names.append("ENSG")

                if 'Ensembl Transcript ID (ENST)' in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                             index_reference_transcript].ENST)
                    column_names.append("Ref_transcript_ID")
                    column_values.append(transcript.ENST)
                    column_names.append("Iso_transcript_ID")

                if 'Ensembl Protein ID (ENSP)' in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[index_reference_transcript].ENSP)
                    column_names.append("Ref_protein_ID")
                    column_values.append(transcript.ENSP)
                    column_names.append("Iso_protein_ID")

                if "Transcript name" in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                             index_reference_transcript].transcript_name)
                    column_names.append("Ref_transcript_name")
                    column_values.append(transcript.transcript_name)
                    column_names.append("Iso_transcript_name")

                if 'Refseq Gene ID (Number)' in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].refseq_gene_ID)
                    column_names.append("Gene_name")

                column_names = column_names + ['AA', 'ReferencePos', 'IsoformPos']
                return column_values, column_names

            #save the results of each alignment in a list outside the loop
            for indexiterator in range(0, len(aminoacids)):
                column_values, column_names = get_selected_columns_attributes_and_column_names(chosen_columns)
                positions = [aminoacids[indexiterator], reference_position_list[indexiterator],isoform_positions_list[indexiterator]]
                nested_list_alignment = column_values + positions
                list_of_all_alignments.append(nested_list_alignment)

        #if multiple ID's list are returned and the pandas dataframe is created later in the create_table_for_dict_of_gene_objects
        if "column_names" in locals(): #column_values, column_names are only generated if len(aminoacid) >0, which means that there has to be at least one match
            if one_ID:
                df = pd.DataFrame(list_of_all_alignments, columns=(column_names))
                return df
            else:
                return list_of_all_alignments, column_names
        else:
            return ('no','matches')


    @staticmethod
    def create_table_for_dict_of_gene_objects(nested_dict, list_of_gene_objects, chosen_columns, match, mismatch,
                                              open_gap_penalty, gap_extension_penalty, exon_length_AA, ID_type="ENSP"):
        list_of_alignments = []
        for gene in nested_dict.items():
            index_of_gene = list(gene[1].keys())[0]
            index_of_reference_transcript = list(gene[1].values())[0]
            if len(list_of_gene_objects[
                       index_of_gene].protein_sequence_isoform_collection) > 1:  # check if there is even more than one isoform
                list_of_dataframe, column_names = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                   list_of_gene_objects, index_of_gene,
                                                                                   chosen_columns, match, mismatch,
                                                                                   open_gap_penalty,
                                                                                   gap_extension_penalty,
                                                                                   exon_length_AA, ID_type=ID_type,
                                                                                   one_ID=False)

                if (list_of_dataframe, column_names) != ('no','matches'): #don't add gene object alignments with no matches at all
                    list_of_alignments = list_of_alignments + list_of_dataframe

        df = pd.DataFrame(list_of_alignments, columns=(column_names))
        return df
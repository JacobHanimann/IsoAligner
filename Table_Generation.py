import pandas as pd
from Alignment import *

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
    def create_table_for_one_gene_object(chosen_reference, list_of_gene_objects, index_of_gene, chosen_columns, match,
                                         mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA,
                                         ID_type="ENSP", one_ID=True):
        '''
        function to create a pandas dataframe for only one gene object
        :param chosen_reference:
        :param list_of_gene_objects:
        :param index_of_gene:
        :param chosen_columns:
        :return: pandas dataframe
        '''

        list_of_all_alignments = []
        if one_ID:  # chosen_reference is a transcript name
            for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
                if getattr(transcript, ID_type) == chosen_reference:
                    reference_protein_sequence = getattr(transcript, "Protein_sequence")
                    break
        else:
            reference_protein_sequence = list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                chosen_reference].Protein_sequence  # chosen_reference is an index

        for index, transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
            if one_ID:
                if getattr(transcript, ID_type) == chosen_reference:
                    index_reference_transcript = index
                    continue
            else:
                if index == chosen_reference:
                    index_reference_transcript = index
                    continue

            aminoacids, reference_position_list, isoform_positions_list = Alignment.map_AA_Needleman_Wunsch_with_exon_check(
                reference_protein_sequence, transcript.Protein_sequence, match, mismatch, open_gap_penalty,
                gap_extension_penalty, exon_length_AA)[1:4]

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

                if "Ensembl Gene ID" in chosen_columns:
                    column_values.append(transcript.ENSG)
                    column_names.append("ENSG")

                if "Ensembl Transcript ID" in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                             index_reference_transcript].ENST)
                    column_names.append("Ref_transcript_ID")
                    column_values.append(transcript.ENST)
                    column_names.append("Iso_transcript_ID")

                if "Ensembl Protein ID" in chosen_columns:
                    column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                             index_reference_transcript].ENSP)
                    column_names.append("Ref_protein_ID")
                    column_values.append(transcript.ENSP)
                    column_names.append("Iso_protein_ID")

                column_names = column_names + ['AA', 'ReferencePos', 'IsoformPos']
                return column_values, column_names

            for indexiterator in range(0, len(aminoacids)):
                column_values, column_names = get_selected_columns_attributes_and_column_names(chosen_columns)
                positions = [aminoacids[indexiterator], reference_position_list[indexiterator],
                             isoform_positions_list[indexiterator]]
                nested_list_alignment = column_values + positions
                list_of_all_alignments.append(nested_list_alignment)

        if one_ID:
            df = pd.DataFrame(list_of_all_alignments, columns=(column_names))
            return df
        else:
            return list_of_all_alignments, column_names



    @staticmethod
    def create_table_for_dict_of_gene_objects(nested_dict, list_of_gene_objects, chosen_columns, match, mismatch,
                                              open_gap_penalty, gap_extension_penalty, exon_length_AA, ID_type="ENSP"):
        list_of_alignments = []
        for gene in nested_dict.items():
            index_of_gene = list(gene[1].keys())[0]
            index_of_reference_transcript = list(gene[1].values())[0]
            if len(list_of_gene_objects[
                       index_of_gene].protein_sequence_isoform_collection) > 1:  # check if there is even more than one isoform
                list_of_dataframe, column_names = create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                   list_of_gene_objects, index_of_gene,
                                                                                   chosen_columns, match, mismatch,
                                                                                   open_gap_penalty,
                                                                                   gap_extension_penalty,
                                                                                   exon_length_AA, ID_type=ID_type,
                                                                                   one_ID=False)
                list_of_alignments = list_of_alignments + list_of_dataframe

        df = pd.DataFrame(list_of_alignments, columns=(column_names))
        return df


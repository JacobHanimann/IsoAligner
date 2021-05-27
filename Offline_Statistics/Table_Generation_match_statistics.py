import pandas as pd
from Alignment import *
import streamlit as st
from Statistics import *
import Bio
from Bio.pairwise2 import format_alignment, align

class Table_Generation_match:
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
    def get_selected_columns_attributes_and_column_names(chosen_columns,index_of_gene, index_reference_transcript,transcript,list_of_gene_objects):
        '''
        function to get the column values selected for each row
        :param chosen_columns: from multiselectbox
        '''
        column_values = []
        column_names = []

        if 'All Columns' in chosen_columns:
            chosen_columns = ['Gene name', 'Ensembl Gene ID (ENSG)', 'Ensembl Transcript ID (ENST)', 'Ensembl Protein ID (ENSP)', 'Transcript name',
             'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)','Refseq Protein ID (NP)','UCSC Stable ID (uc)','Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID', 'Uniparc ID',
             'Ensembl Gene ID version (ENSG.Number)', 'Ensembl Transcript ID version (ENST.Number)', 'Ensembl Protein ID version (ENSP.Number)', 'Refseq Transcript ID version (NM.Number)', 'Refseq Transcript ID version (NP.Number)',
             'HGNC ID (HGNC:Number)']

        if "Gene name" in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].ensembl_gene_symbol)
            column_names.append("Gene_name")

        if 'HGNC ID (HGNC:Number)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].HGNC)
            column_names.append("HGNC_ID")

        if 'Refseq Gene ID (Number)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].refseq_gene_ID)
            column_names.append("Refseq_Gene_ID")

        if "Ensembl Gene ID (ENSG)" in chosen_columns:
            column_values.append(transcript.ENSG)
            column_names.append("ENSG")
        if 'Ensembl Gene ID version (ENSG.Number)' in chosen_columns:
            column_values.append(transcript.ENSG_version)
            column_names.append("ENSG_version")

        if 'Ensembl Transcript ID (ENST)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].ENST)
            column_names.append("Ref_transcript_ID")
            column_values.append(transcript.ENST)
            column_names.append("Iso_transcript_ID")
        if 'Ensembl Transcript ID version (ENST.Number)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].ENST_version)
            column_names.append("Ref_ts_ID_version")
            column_values.append(transcript.ENST_version)
            column_names.append("Iso_ts_ID_version")

        if 'Ensembl Protein ID (ENSP)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].ENSP)
            column_names.append("Ref_protein_ID")
            column_values.append(transcript.ENSP)
            column_names.append("Iso_protein_ID")
        if 'Ensembl Protein ID version (ENSP.number)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].ENSP_version)
            column_names.append("Ref_prot_ID_ver")
            column_values.append(transcript.ENSP_version)
            column_names.append("Iso_prot_ID_ver")

        if "Transcript Name" in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].transcript_name)
            column_names.append("Ref_transcript_name")
            column_values.append(transcript.transcript_name)
            column_names.append("Iso_transcript_name")

        if 'Refseq Transcript ID (NM)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].refseq_NM)
            column_names.append("Ref_NM_ID")
            column_values.append(transcript.refseq_NM)
            column_names.append("Iso_NM_ID")
        if 'Refseq Transcript ID version (NM.Number)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].refseq_NM_version)
            column_names.append("Ref_NM_ID_ver")
            column_values.append(transcript.refseq_NM_version)
            column_names.append("Iso_NM_ID_ver")

        if 'Refseq Protein ID (NP)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].refseq_NP)
            column_names.append("Ref_NP_ID")
            column_values.append(transcript.refseq_NP)
            column_names.append("Iso_NP_ID")
        if 'Refseq Protein ID version (NP.Number)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].refseq_NP_version)
            column_names.append("Ref_NP_ID")
            column_values.append(transcript.refseq_NP_version)
            column_names.append("Iso_NP_ID")

        if 'UCSC Stable ID (uc)' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].UCSC_stable_ID)
            column_names.append("Ref_UCSC_ID")
            column_values.append(transcript.UCSC_stable_ID)
            column_names.append("Iso_UCSC_ID")

        if 'Uniprot Name ID' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].uniprot_name_ID)
            column_names.append("Uniprot_name_ID")

        if 'Uniprot Accession ID' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].uniprot_accession)
            column_names.append("Ref_uniprot_accession")
            column_values.append(transcript.uniprot_accession)
            column_names.append("Iso_uniprot_accession")
        if 'Uniparc ID' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].uniprot_uniparc)
            column_names.append("Ref_uniprot_uniparc")
            column_values.append(transcript.uniprot_uniparc)
            column_names.append("Iso_uniprot_uniparc")
        if 'Uniprot Isoform ID' in chosen_columns:
            column_values.append(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
                                     index_reference_transcript].uniprot_isoform)
            column_names.append("Ref_uniprot_isoform")
            column_values.append(transcript.uniprot_isoform)
            column_names.append("Iso_uniprot_isoform")

        column_names = column_names + ['AA', 'ReferencePos', 'IsoformPos']
        return column_values, column_names

    @staticmethod
    def create_table_for_one_gene_object(index_reference_transcript, list_of_gene_objects, index_of_gene, chosen_columns, match,
                                         mismatch, open_gap_penalty, gap_extension_penalty,
                                         one_ID=True,exon_length=None):
        '''
        one_ID: function returns pandas dataframe directly of alignments
        multiple IDs: function returns with list of the computed alignments which is then further used in create_table_for_dict_of_gene_objects function

        :param chosen_reference: index
        :param list_of_gene_objects: (library)
        :param index_of_gene:
        :param chosen_columns: (chosen from selectbox)
        :return: pandas dataframe or lists of computed alignment
        '''
        list_of_all_alignments = []
        gene_check_list = []

        #select reference protein sequence
        reference_protein_sequence = list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[index_reference_transcript].protein_sequence  # chosen_reference is an index (table generation for multiple ID's)

        #set the minimal exon length
        #by the user when looking at one gene
        if exon_length !=None and one_ID:
            exon_length_AA = exon_length
        #checking library for value
        elif list_of_gene_objects[index_of_gene].minimal_exon_length!= None:
            if list_of_gene_objects[index_of_gene].minimal_exon_length >4:
                exon_length_AA = list_of_gene_objects[index_of_gene].minimal_exon_length
            else:
                exon_length_AA = 11
        else:
            #if there is no value in the library
            exon_length_AA = 11

            #create alignment for each alternative_ID isoform
        aa_correct = 0
        aa_false = 0
        for index, transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
            if index == index_reference_transcript: #do not align the reference transcript with itself
                continue
            aminoacids, reference_position_list, isoform_positions_list,isoform_pattern_check = Alignment.map_AA_Needleman_Wunsch_with_exon_check( #add isoform_check_list and [1:5] and uncomment all lines with aa_correct, aa_false associations for false,positive statistics
                reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)[1:5]

            correct, false = Statistics.isoform_form_check_stats(isoform_pattern_check)
            aa_correct = aa_correct + correct
            aa_false = aa_false + false


        return (aa_correct, aa_false)


    @staticmethod
    #@st.cache()
    #@st.cache(hash_funcs={Alignment.map_AA_Needleman_Wunsch_with_exon_check: hash})
    def create_table_for_dict_of_gene_objects(nested_dict, list_of_gene_objects, chosen_columns, match, mismatch,
                                              open_gap_penalty, gap_extension_penalty, conventional=False):

        if conventional:
            open_gap_penalty=0
        list_of_alignments = []
        total_genes = len(nested_dict)
        correct_aa = 0
        false_aa = 0
        with st.spinner('Generating Mapping Table . . . '):
            my_bar = st.progress(0.0)
            for count,gene in enumerate(nested_dict.items()):
                percent = int(round(100*(count+1)/total_genes,1))
                print(percent,"%")
                index_of_gene = list(gene[1].keys())[0]
                index_of_reference_transcript = list(gene[1].values())[0]
                if len(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection) > 1:  # check if there is even more than one isoform
                    gene_isoform_check_list = Table_Generation_match.create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                       list_of_gene_objects, index_of_gene,
                                                                                       chosen_columns, match, mismatch,
                                                                                       open_gap_penalty,
                                                                                       gap_extension_penalty,
                                                                                       one_ID=False)

                    if gene_isoform_check_list != ('no','matches'): #don't add gene object alignments with no matches at all
                        #list_of_alignments = list_of_alignments + list_of_dataframe
                        #my_bar.progress(percent)
                        correct_aa = correct_aa + gene_isoform_check_list[0]
                        false_aa = false_aa + gene_isoform_check_list[1]
                    else:
                        print(gene)
            #my_bar.empty()
            #df = pd.DataFrame(list_of_alignments, columns=(column_names))

        return correct_aa, false_aa

    @staticmethod
    def display_filter_option_AA():
        option,blankspace = st.beta_columns([1,2])
        with option:
            value = st.text_input('Filter mapping table for specific value', value='')
        return value

    @staticmethod
    def filter_all_columns_of_df(value, dataframe):
        columns = dataframe.columns
        if value.isnumeric():
            value = int(value)
        final_dataframe = pd.DataFrame(columns=columns)
        for column in columns:
            filter_df = dataframe.loc[(dataframe[column] == value)]
            if not filter_df.empty:
                final_dataframe = pd.concat([final_dataframe, filter_df], ignore_index=True).drop_duplicates()
            else:
                pass
        return final_dataframe
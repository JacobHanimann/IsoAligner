import pandas as pd
from IsoAligner_core.Alignment import *
from Streamlit_app.Statistics import *


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
                exon_length_AA = 10
        else:
            #if there is no value in the library
            exon_length_AA = 10

            #create alignment for each alternative_ID isoform
        for index, transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
            if index == index_reference_transcript: #do not align the reference transcript with itself
                continue
            aminoacids, reference_position_list, isoform_positions_list = Alignment.map_AA_Needleman_Wunsch_with_exon_check( #add isoform_check_list and [1:5] and uncomment all lines with aa_correct, aa_false associations for false,positive statistics
                reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)[1:4]


            #save the results of each alignment in a list outside the loop
            for indexiterator in range(0, len(aminoacids)):
                column_values, column_names = Table_Generation.get_selected_columns_attributes_and_column_names(chosen_columns,index_of_gene, index_reference_transcript,transcript,list_of_gene_objects)
                positions = [aminoacids[indexiterator], reference_position_list[indexiterator],isoform_positions_list[indexiterator]]
                nested_list_alignment = column_values + positions
                list_of_all_alignments.append(nested_list_alignment)


        #if multiple ID's list are returned and the pandas dataframe is created later in the create_table_for_dict_of_gene_objects
        if "column_names" in locals(): #column_values, column_names are only generated if len(aminoacid) >0, which means that there has to be at least one match
            if one_ID:
                df = pd.DataFrame(list_of_all_alignments, columns=(column_names))
                return df

            else:
                return list_of_all_alignments, column_names, #(aa_correct, aa_false)
        else:
            return ('no','matches')


    @staticmethod
    #@st.cache()
    #@st.cache(hash_funcs={Alignment.map_AA_Needleman_Wunsch_with_exon_check: hash})
    def create_table_for_dict_of_gene_objects(nested_dict, list_of_gene_objects, chosen_columns, match, mismatch,
                                              open_gap_penalty, gap_extension_penalty):

        list_of_alignments = []
        total_genes = len(nested_dict)
        #correct_aa = 0
        #false_aa = 0
        with st.spinner('Generating Mapping Table . . . '):
            my_bar = st.progress(0.0)
            for count,gene in enumerate(nested_dict.items()):
                percent = int(round(100*(count+1)/total_genes,1))
                print(percent)
                index_of_gene = list(gene[1].keys())[0]
                index_of_reference_transcript = list(gene[1].values())[0]
                if len(list_of_gene_objects[
                           index_of_gene].protein_sequence_isoform_collection) > 1:  # check if there is even more than one isoform
                    list_of_dataframe, column_names = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                       list_of_gene_objects, index_of_gene,
                                                                                       chosen_columns, match, mismatch,
                                                                                       open_gap_penalty,
                                                                                       gap_extension_penalty,
                                                                                       one_ID=False)

                    if (list_of_dataframe, column_names) != ('no','matches'): #don't add gene object alignments with no matches at all
                        list_of_alignments = list_of_alignments + list_of_dataframe
                        my_bar.progress(percent)

                    else:
                        print('isoforms with no overlaps apparently:', gene)
            my_bar.empty()
            if list_of_alignments:
                df = pd.DataFrame(list_of_alignments, columns=(column_names))
            else:
                return ('no','matches')

        return df

    @staticmethod
    def display_filter_option_AA():
        option,blankspace = st.columns([1,2])
        with option:
            value = st.text_input('Filter mapping table for specific value:', value='')
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
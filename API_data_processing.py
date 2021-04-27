from flask import Flask
from flask_restful import Api, Resource, reqparse, abort, fields, marshal_with
from flask_caching import Cache
from Visualise_Alignment import *
from Alignment import *
import pickle
from Gene import *
from Protein_isoform import *
import sys
# insert at position 1 in the path, as 0 is the path of this file.
#sys.path.insert(1, '../')
import sys
from Gene import *
from Protein_isoform import *
from Streamlit_community import *
from Input_flow import *
from Streamlit_Pop_ups import *
from Alignment import *
from Visualise_Alignment import *
from User_Input_Preparation import *
from Input_flow import *
from Table_Generation import *
from PIL import Image
from Statistics import *
from Table_Generation import *
from flask_caching import Cache

class Data_processing():
    pass


    @staticmethod
    def align_sequences(input1, input2):
        needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(input1, input2, 1, -2, -1.75, 0, 5)
        isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
        percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check, input1, input2)
        alignment_string = Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta,alignment_isoform_fasta,isoform_pattern_check,percentage_reference, percentage_isoform)
        return alignment_string

    @staticmethod
    def search_and_generate_nested_dict(ID_to_search, list_of_gene_objects):
        dict_of_IDs = Input_preparation.identify_IDs_from_user_text_input(ID_to_search)
        gene_index = Input_flow.search_through_database_with_known_ID_Type(list_of_gene_objects, dict_of_IDs)
        verified_gene_index = Input_flow.remove_dict_elements_with_no_gene_object_match(gene_index)
        if not verified_gene_index:
            return None
        nested_dict = Input_flow.generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs,verified_gene_index,list_of_gene_objects)
        return nested_dict,dict_of_IDs


    @staticmethod
    def choose_column_if_two_IDs(reference_ID, alternative_ID):
        attribute_column_dict = {'gene_name':'Gene name','ENSG':'Ensembl Gene ID (ENSG)', 'ENST':'Ensembl Transcript ID (ENST)',
                      'ENSP':'Ensembl Protein ID (ENSP)', 'transcript_name':'Transcript name',
                      'refseq_gene_ID':'Refseq Gene ID (Number)', 'refseq_NM':'Refseq Transcript ID (NM)','refseq_NP': 'Refseq Protein ID (NP)',
                      'UCSC_stable_ID':'UCSC Stable ID (uc)',
                      'uniprot_name_ID':'Uniprot Name ID','uniprot_accession':'Uniprot Accession ID','uniprot_isoform': 'Uniprot Isoform ID','uniprot_uniparc': 'Uniparc ID',
                      'ENSG_version':'Ensembl Gene ID version (ENSG.Number)','ENST_version':'Ensembl Transcript ID version (ENST.Number)',
                      'ENSP_version':'Ensembl Protein ID version (ENSP.Number)', 'refseq_NM_version':'Refseq Transcript ID version (NM.Number)',
                      'refseq_NP_version':'Refseq Transcript ID version (NP.Number)',
                      'HGNC':'HGNC ID (HGNC:Number)'}
        chosen_column = [attribute_column_dict[reference_ID]+attribute_column_dict[alternative_ID]]
        return chosen_column



    @staticmethod
    def choose_mapping_table_columns(table_ids,reference_id, alternative_id,two_IDs=False):
        chosen_columns = ['Gene name']
        if table_ids!=None:
            if 'ensembl' in table_ids:
                chosen_columns = chosen_columns + ['Ensembl Gene ID (ENSG)', 'Ensembl Transcript ID (ENST)',
                          'Ensembl Protein ID (ENSP)', 'Transcript name', 'Ensembl Gene ID version (ENSG.Number)', 'Ensembl Transcript ID version (ENST.Number)']
            if 'refseq' in table_ids:
                chosen_columns = chosen_columns + ['Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)', 'Refseq Protein ID (NP)','Refseq Transcript ID version (NM.Number)',
                          'Refseq Transcript ID version (NP.Number)']
            if 'uniprot' in table_ids:
                chosen_columns = chosen_columns + ['Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID', 'Uniparc ID']
            if 'ucsc' in table_ids:
                chosen_columns = chosen_columns + [ 'UCSC Stable ID (uc)']
            if 'hgnc' in table_ids:
                chosen_columns = chosen_columns + ['HGNC ID (HGNC:Number)']

        if table_ids==None or len(chosen_columns)==1 or 'all' in table_ids:
            chosen_columns = ['Gene name', 'Ensembl Gene ID (ENSG)', 'Ensembl Transcript ID (ENST)',
                      'Ensembl Protein ID (ENSP)', 'Transcript name',
                      'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)', 'Refseq Protein ID (NP)',
                      'UCSC Stable ID (uc)',
                      'Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID', 'Uniparc ID',
                      'Ensembl Gene ID version (ENSG.Number)', 'Ensembl Transcript ID version (ENST.Number)',
                      'Ensembl Protein ID version (ENSP.Number)', 'Refseq Transcript ID version (NM.Number)',
                      'Refseq Transcript ID version (NP.Number)',
                      'HGNC ID (HGNC:Number)']
        if two_IDs==True and table_ids==None:
            chosen_columns = Data_processing.choose_column_if_two_IDs(reference_id,alternative_id)
        return chosen_columns


    @staticmethod
    def create_mapping_table_of_two_IDs(list_of_gene_objects,index_of_gene,index_reference_transcript,index_alternative_transcript,chosen_columns, match, mismatch,open_gap_penalty, gap_extension_penalty,exon_length_AA):
        list_of_all_alignments = []

        # reference protein sequence
        reference_protein_sequence = list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
            index_reference_transcript].protein_sequence

        # alternative protein sequence
        alternative_protein_sequence = list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[
            index_alternative_transcript].protein_sequence

        #protein isoform object of alternative sequence
        transcript = list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[index_alternative_transcript]

        aminoacids, reference_position_list, isoform_positions_list = Alignment.map_AA_Needleman_Wunsch_with_exon_check(reference_protein_sequence, alternative_protein_sequence, match, mismatch, open_gap_penalty,
                                                                                            gap_extension_penalty, exon_length_AA)[1:4]

        # save the results of each alignment in a list outside the loop
        for indexiterator in range(0, len(aminoacids)):
            column_values, column_names = Table_Generation.get_selected_columns_attributes_and_column_names(
                chosen_columns, index_of_gene, index_reference_transcript, transcript, list_of_gene_objects)
            positions = [aminoacids[indexiterator], reference_position_list[indexiterator],
                         isoform_positions_list[indexiterator]]
            nested_list_alignment = column_values + positions
            list_of_all_alignments.append(nested_list_alignment)

        df = pd.DataFrame(list_of_all_alignments, columns=(column_names))
        return df


    @staticmethod
    def extract_specific_position_mapping_table(df,AA_position):
        '''
        function that extracts specified aminoacid postion from one to another isoform
        :param df: pandas dataframe
        :param position: aminoacid symbol and position, for example: A322
        :return: amino acid position of the alternative isoform
        '''
        match_list = re.findall('\d+',AA_position)
        if not match_list:
            return 'aa_position must be an integer'
        position_number = int(AA_position)
        #df.loc[(df['AA'] == AA) & (df['ReferencePos'] == position_number)]
        row = df.loc[df['ReferencePos'] == position_number]
        if not row.empty:
            AA_position_new = str(list(row['IsoformPos'])[0])
            return AA_position_new
        else:
            return 'position not found on reference isoform', 400
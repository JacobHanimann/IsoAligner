import pandas as pd
from IsoAligner_core.Alignment import *
from Streamlit_app.Statistics import *

class Table_Generation_match:
    pass


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
            if list_of_gene_objects[index_of_gene].minimal_exon_length >=3:
                exon_length_AA = list_of_gene_objects[index_of_gene].minimal_exon_length
            else:
                exon_length_AA = 10
        else:
            #if there is no value in the library
            exon_length_AA = 10

            #create alignment for each alternative_ID isoform
        aa_correct = 0
        aa_false = 0
        aa_mismatch =0
        for index, transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
            if index == index_reference_transcript: #do not align the reference transcript with itself
                continue
            aminoacids, reference_position_list, isoform_positions_list,isoform_pattern_check = Alignment.map_AA_Needleman_Wunsch_with_exon_check( #add isoform_check_list and [1:5] and uncomment all lines with aa_correct, aa_false associations for false,positive statistics
                reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)[1:5]

            correct, false,mmatch = Statistics.isoform_form_check_stats(isoform_pattern_check)
            aa_correct = aa_correct + correct
            aa_false = aa_false + false
            aa_mismatch = aa_mismatch +mmatch


        return (aa_correct, aa_false, aa_mismatch)



    @staticmethod
    def create_table_for_dict_of_gene_objects(nested_dict, list_of_gene_objects, chosen_columns, match, mismatch,
                                              open_gap_penalty, gap_extension_penalty, conventional=False):

        if conventional:
            open_gap_penalty=0
        list_of_alignments = []
        total_genes = len(nested_dict)
        correct_aa = 0
        false_aa = 0
        mismatch_aa = 0
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
                        mismatch_aa = mismatch_aa +gene_isoform_check_list[2]
                    else:
                        print(gene)
            #my_bar.empty()
            #df = pd.DataFrame(list_of_alignments, columns=(column_names))

        return correct_aa, false_aa, mismatch_aa
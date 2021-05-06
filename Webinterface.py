import streamlit as st
import sys
import SessionState
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


#declare session state variables
ss = SessionState.get(clicked=False,searched_clicked=False, align_clicked=False, generate=False,run_id=0,example=False, clear_button=False)

#import database
list_of_gene_objects = Input_flow.import_data_from_github('list_of_gene_objects_4th_may.txt.gz')


#Streamlit website
def main():
    """ Isoform Alignment Tool """
    #Background
    #Streamlit_community.set_png_as_page_bg('TransferMessengerRNAt.png')

    #Sidebar
    activity = ['Alignment Tool', 'REST API & Downloads', 'Manual & About']
    st.sidebar.markdown("## Navigation")
    choice = st.sidebar.radio("Go to", activity)
    st.sidebar.write("\n")

    #Alignment tool section
    if choice == 'Alignment Tool':
        header, tRNA = st.beta_columns([3, 1.1])
        with header:
            # Title
            st.title(" Amino Acid Isoform Aligner")
            st.subheader("Map Amino Acid Positions Across Isoforms")
            st.write(
                "Align isoforms dynamically with the Needleman-Wunsch algorithm and set the minimal exon length to discard falsely mapped positions."
                " The current human isoform library consists of ~18k protein coding genes covering ~130k protein sequences and ~1.3M mapped isoform ID's from Ensembl, Uniprot, Refseq and HGNC.")
        with tRNA:
            st.write('\n')
            #st.markdown('[this is a text link](upload://7FxfXwDqJIZdYJ2QYADywvNRjB.png)')
            #st.markdown('[![this is an image link](upload://TransferMessengerRNA.tif)](https://streamlit.io)')
            image2 = Image.open('Pictures/Spliceosome_yeast_small.tif')
            st.image(image2,use_column_width=True)

        st.write("--------------------------")
        st.sidebar.markdown("### 🧬️Organism")
        st.sidebar.selectbox('Select species', ['🧍🏽Homo Sapiens', '🐁 Mouse (next release)'])
        st.write('\n')

        #fixed in put area
        title, example_button = st.beta_columns([3.5,1])
        with title:
            st.markdown("#### Input")
        with example_button:
            if st.button('Show Example'):
                ss.example = True
                ss.searched_clicked = False


        if ss.example:
            #input1 = st.text_area('Paste multiple ID\'s comma or newline and click on search library for ID\'s. Go to "Manual" for further information', '''EGFR, Q9Y6I3, Q9Y6I3-1, ENSG00000074410, ENSG00000164690.2, ENSP00000005756, ENSP00000075430.7, HGNC:10728, UPI00022F85F1, uc060zgm.1, NP_004702.2, ENST00000554846.5, SEMA5B-211, BRAF_HUMAN, NM_001304833, 1232, ENST00000551640, NP_775733, NM_003769.3''',key=ss.run_id)
            input1 = st.text_area('Paste multiple ID\'s comma or newline and click on search library for ID\'s. Go to "Manual" for further information',
                                  '''EGFR, Q9Y6I3 - 1, ENSG00000074410, ENSP00000075430.7, HGNC: 10728, UPI00022F85F1, NP_004702.2, SEMA5B - 211, BRAF_HUMAN, NM_001304833''',key=ss.run_id)
        else:
            input1 = st.text_area('Paste any Ensembl/Uniprot/Refseq ID\'s, gene names or a raw amino acid sequence: ', '''''',key=ss.run_id)
        file_upload, search_button = st.beta_columns([2,1])
        with file_upload:
            file_wanted = st.checkbox("upload list of ID's or gene names")
            if file_wanted:
                input1 = st.file_uploader("Accepted ID's: Ensembl, Refseq, HGNC, Uniprot (Accession/Isoform/Uniparc)", type=[ "txt"])
                if input1 is not None:
                    input1 = input1.getvalue().decode("UTF-8")
                    ss.searched_clicked = True
            raw_aa = st.checkbox("insert 2nd raw amino acid manually")
        with search_button:
            search = st.button('Search Library for ID\'s')
            if search:
                ss.searched_clicked =True

        #set default for displaying second text_area input for input2
        using_IDs= False

        if ss.searched_clicked:
            dict_of_IDs = Input_preparation.identify_IDs_from_user_text_input(input1)
            #st.write(dict_of_IDs)
            #if an element is an amino acid:
                #warning message
            #else (no AA sequence)
            input1_IDs = Input_flow.search_through_database_with_known_ID_Type(list_of_gene_objects, dict_of_IDs)
            #st.write(input1_IDs)
            Input_flow.show_which_elements_were_not_found(input1_IDs)
            cleaned_input1_IDs=Input_flow.remove_dict_elements_with_no_gene_object_match(input1_IDs)
            #execute nested dict only if dict is still existent...
            nested_dict = Input_flow.generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs, cleaned_input1_IDs,list_of_gene_objects)
            #st.write(nested_dict)
            no_elements = False
            if len (nested_dict)==0:
                no_elements = True

        #case of using one ID
        if ss.searched_clicked and not no_elements and  len(input1_IDs) == 1: #check if dictionary is not empty
            using_IDs = True
            #st.write(input1_IDs)
            #st.write(list(input1_IDs.values())[0])
            st.markdown("### Alignments")
            reference_select, number_of_entries = st.beta_columns([1,2.5])
            with reference_select:
                chosen_gene = list(nested_dict.keys())[0]
                index_gene_object = list(list(nested_dict.values())[0].keys())[0]
                transcript_list,index_gene = Visualise_Alignment.fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,nested_dict,chosen_gene)
                chosen_reference = st.selectbox('Choose your reference transcript: ',[transcript[0] for transcript in transcript_list])
                index_of_reference_transcript = Visualise_Alignment.get_index_of_chosen_transcript(chosen_reference,transcript_list)
            with number_of_entries:
                st.write('\n')
                st.write('\n')
                st.write('\n')
                st.write('Number of entries:',len(transcript_list))
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, index_gene_object)
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
            st.write('\n')
            st.text('\n')
            Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript,list_of_gene_objects,index_gene_object,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            #Table section
            chosen_columns = Input_flow.chose_columns()
            generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,list_of_gene_objects,index_gene_object,chosen_columns,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length=exon_length_AA)
            st.text('\n')
            st.write(generated_table)
            st.text('\n')
            Input_flow.generate_download_section(generated_table)


        #case of using multiple ID's
        elif ss.searched_clicked and len(input1_IDs) > 1 and not no_elements:
            using_IDs = True
            st.markdown("### Alignments")
            st.text('\n')
            genes, reference = st.beta_columns([2,1.7])
            with genes:
                chosen_gene = st.selectbox('Select Gene',Visualise_Alignment.create_list_gene_selection(list_of_gene_objects,nested_dict))
            with reference:
                transcript_list, index_gene = Visualise_Alignment.fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,nested_dict, chosen_gene)
                chosen_reference = st.selectbox('Choose your reference transcript: ',[transcript[0] for transcript in transcript_list])
                index_of_reference_transcript = Visualise_Alignment.get_index_of_chosen_transcript(chosen_reference,transcript_list)
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, index_gene)
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
            st.write('\n')
            st.text('\n')
            gene_index = list(nested_dict[re.split(' \(',Visualise_Alignment.clean_chosen_gene(chosen_gene))[0]])[0]
            #st.write('indexes of gene objects:')
            #st.write(nested_dict)
            Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript, list_of_gene_objects,gene_index, match, mismatch, open_gap_penalty, gap_extension_penalty,exon_length_AA)
            # Table section
            chosen_columns = Input_flow.chose_columns()
            df_all = Table_Generation.create_table_for_dict_of_gene_objects(nested_dict,list_of_gene_objects,chosen_columns, match, mismatch, open_gap_penalty, gap_extension_penalty)
            if not df_all.empty:
                st.write(df_all)
                st.text('\n')
                Input_flow.generate_download_section(df_all)
            else:
                st.warning('No amino acid positions mapped currently. Tweak function parameters to generate matches.')


        #Input 2 Area
        if using_IDs== False and raw_aa:
            input2 = st.text_area('Paste Amino Acid sequence of alternative_ID isoform: ', '''''', key=ss.run_id)
            align=st.button('Align')
            if align:
                ss.align_clicked = True
                ss.searched_clicked = False
            if input1 != "" and input2 != "" and ss.align_clicked and ss.searched_clicked==False:
                #Sidebar pop up, make function out of it?
                match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, index_gene, raw=True)
                st.markdown("### Results")
                #st.write("\n")
                #st.markdown("##### Unfiltered Alignment:")
                #st.write("\n")
                needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(input1,input2, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)
                isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
                #st.text(Alignment_preview)
                st.write("\n")
                st.markdown("##### Alignment")
                st.write("\n")
                percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check,input1,input2)
                st.text(Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta,alignment_isoform_fasta,isoform_pattern_check,percentage_reference,percentage_isoform))
                st.write("\n")
                st.markdown(" ###### ℹ️Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
                st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
                st.write("\n")
                st.write("\n")
                generated_table = Table_Generation.create_pandas_dataframe_raw_aa_sequence(needleman_mapped)
                if not generated_table.empty:
                    table, whitespace, download = st.beta_columns([1,0.2,1])
                    st.write("\n")
                    with download:
                        st.write("\n")
                        st.write("\n")
                        st.markdown("#### 📁 Download")
                        sep_choice = st.radio('Choose file format:', ['tsv', 'csv'])
                        if sep_choice == "tsv":
                            sep = '\t'
                        else:
                            sep = ','
                        st.markdown(Streamlit_community.get_table_download_link(generated_table, 'dataframe.' + sep_choice, sep),unsafe_allow_html=True)
                    with table:
                        st.markdown("##### Correctly mapped AA positions")
                        st.write("\n")
                        st.write(generated_table)
                else:
                    st.warning('No amino acid positions mapped currently. Tweak function parameters to generate matches.')

        #Clear all button
        st.write("--------------------------")
        placehold, clear_all = st.beta_columns([6.85, 1])
        with clear_all:
           if st.button('Clear All'):
              ss.clear_button = True
              ss.run_id +=1
              ss.example= False
              ss.searched_clicked= False
              Streamlit_community.rerun_script_from_top()



    elif choice == 'REST API & Downloads':
        st.title(" Amino Acid Isoform Aligner")
        st.header("REST API")
        st.write('soon available')
        st.write("--------------------------")
        st.header("Downloads")
        st.write('soon available')
        st.write("--------------------------")

    if choice == 'Manual & About':

        st.title(" Amino Acid Isoform Aligner")
        st.markdown("### Biologically Appropriate Alignment of Isoforms:")
        st.write('The challenge of aligning protein sequences of two isoforms is essentially matching the exons correctly as can be seen in the scheme below.'
                 ' The IsoAligner algorithm exploits the biological characteristics of isoforms by simply confining the parameters of the Needleman Wunsch global alignment to assure positional mapping of biologically corresponding amino acid positions only.  ')
        problem_schema = Image.open('Pictures/Mapping_problem_schema.png')
        st.image(problem_schema,use_column_width=True)
        st.write('To avoid and discard falsely mapped positions of distinct exons (like Exon4 and Exon5) the parameters of the alignment are tweaked as follows:')
        needleman, minimal = st.beta_columns([1, 1.])
        with needleman:
            st.markdown("#### Needleman-Wunsch global alignment")
            st.write(" - Big open gap penalty (Default -2) \n"
                 "- Small extend gap penalty (Default 0)")
        with minimal:
            st.markdown("#### Discard falsely matched positions")
            st.write('- By definition of a minimal exon length (Default 5 AA)')
        st.markdown("### Alignment example:")
        st.write("In general, alignment solutions matching identical exons are preferred. The generated matches are then additionally inspected and verified. Alignment sections only containing partial diffuse mapping are being recognised as random matches and are marked as 'x'.")
        example = Image.open('Pictures/example_schema.png')
        st.image(example, use_column_width=True)
        st.write("--------------------------")
        st.markdown("### Manual Alignment Tool")
        st.write("Quick Start: Click on 'Show Example' and then 'Search Library for IDs' to get a overview.")
        st.write("1. Paste ID's, gene names or raw amino acid sequences"
                 "\n    - The current human library consists of ~24k protein coding genes covering ~190k protein sequences and ~1.5M mapped Isoforms ID's from Ensembl, Uniprot, Refseq & HGNC. Included are Gene, Transcript & Protein ID's of various types."
                 "\n    - Press 'Search Library for IDs' or 'Align' to compute alignments"
                 "\n 2. Tweak function parameters in the sidebar and inspect the alignment previews"
                 "\n    - Set the mininmal exon length (in AA)"
                 "\n""    - Set Needleman-Wunsch parameters"
                "\n 3. Select ID's to be included in the mapping table."
                "\n    - Download dataframe: tab or comma separated (tsv/csv) "
                 )
        st.markdown("#### ⚠️ Important:")
        st.write("As soon as multiple ID's are entered, the reference transcript for the generation of the mapping table is automatically chosen, unless a specific transcript or protein ID is used in the input field."
                 " If Gene names or ID's are given as the input, the isoform with the longest sequence is used as the reference to align against.")
        st.write("--------------------------")
        st.markdown("### Human Isoform Library Overview:")
        st.write('\n')
        total_number_of_genes, total_number_of_isoforms, genes_without_isoforms,Ids_in_total = Statistics.list_of_gene_objects_statistics(list_of_gene_objects)
        st.write('Genes: ',total_number_of_genes)
        st.write('Isoforms in total: ',total_number_of_isoforms)
        st.write('Average number of isoform per gene:', round(total_number_of_isoforms/total_number_of_genes,1))
        st.write("ID's in total:",Ids_in_total)
        st.write("Average number of ID's per protein sequence:", round(Ids_in_total/total_number_of_isoforms,1))
        st.write("Gene object attributes:", Gene.list_of_attributes())
        st.write("Collection of Isoform ID's includes:",Protein_isoform.list_of_attributes())
        #st.write('gene objects without isoforms: ',genes_without_isoforms)
        duplicates_number, genes_without_duplicates, redundant_sequences, genes_with_more_than_one_duplicate = Statistics.check_if_there_are_AA_seq_duplicates(list_of_gene_objects)
        #st.write('Genes with AA seq duplicates: ', duplicates_number)
        #st.write('Genes without AA seq duplicates: ', genes_without_duplicates)
        #st.write('Redundant AA sequences:', redundant_sequences)
        #st.write('Genes with more than one AA seq duplicate: ', genes_with_more_than_one_duplicate)
        #st.write('Ensembl IDs: ')
        #st.write('Refseq IDs: ')
        #st.write('Uniprot IDs: ')
        st.write("--------------------------")
        st.markdown("#### Contact:")
        st.write("\n")
        st.write("Please get in touch for suggestions or to report bugs :)")
        st.text('''Gian Jacob Hanimann\nE-mail: GianJacob.Hanimann@usz.ch\nPhone: +41765596015''')
        st.write('Bioinformatics group: https://clinicalcompbio.org/')

        #st.markdown("#### Functions:")
        #code = '''
        #def transform_uploaded_data_type_accordingly(file):
        #    'uploaded files can be different types of files. A transformation is needed to interpret the data correctly
        #    Type of input: FASTA, FA and TXT
        #    Output type: depends on the case'
        #
        #'''
        #st.code(code, language='python')


#Execution

#Default Needleman- Wunsch Parameters:
match=2
mismatch= -1.75
open_gap_penalty= -1
gap_extension_penalty= 0

#Execution
if __name__ == '__main__':
    main()

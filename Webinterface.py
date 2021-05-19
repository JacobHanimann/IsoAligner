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
import time


#declare session state variables
ss = SessionState.get(clicked=False,searched_clicked=False, align_clicked=False, generate=False,run_id=0,example=False, clear_button=False,run_id_table=1, parameters=[1,2,3,4,5])

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


    #Alignment tool section
    if choice == 'Alignment Tool':
        header, tRNA = st.beta_columns([3, 1.1])
        with header:
            # Title
            st.title(" Amino Acid Isoform Aligner")
            st.subheader("Map Amino Acid Positions Across Isoforms")
            st.write(
                "Align protein isoforms interactively with a fitted Needleman-Wunsch algorithm and set the minimal exon length to discard falsely mapped positions."
                " The current human isoform library consists of ~18K protein coding genes covering ~130K protein sequences and +1.3M mapped isoform ID's from Ensembl, Uniprot, Refseq, UCSC and HGNC.")
        with tRNA:
            st.write('\n')
            #st.markdown('[this is a text link](upload://7FxfXwDqJIZdYJ2QYADywvNRjB.png)')
            #st.markdown('[![this is an image link](upload://TransferMessengerRNA.tif)](https://streamlit.io)')
            image2 = Image.open('Pictures/Spliceosome_yeast_small.tif')
            st.image(image2,use_column_width=True)

        st.sidebar.markdown("### 🧬️Organism")
        st.sidebar.selectbox('Select species', ['🧍🏽Homo Sapiens', '🐁 Mouse (next release)'])

        #fixed in put area
        title, example_button = st.beta_columns([3.85,1])
        with title:
            st.markdown("### Input")
        with example_button:
            if st.button('Show Example'):
                ss.example = True
                ss.searched_clicked = False


        if ss.example:
            #input1 = st.text_area('Paste multiple ID\'s comma or newline and click on search library for ID\'s. Go to "Manual" for further information', '''EGFR, Q9Y6I3, Q9Y6I3-1, ENSG00000074410, ENSG00000164690.2, ENSP00000005756, ENSP00000075430.7, HGNC:10728, UPI00022F85F1, uc060zgm.1, NP_004702.2, ENST00000554846.5, SEMA5B-211, BRAF_HUMAN, NM_001304833, 1232, ENST00000551640, NP_775733, NM_003769.3''',key=ss.run_id)
            input1 = st.text_area('Paste multiple ID\'s (comma or newline separated) and click on search library. Go to "Manual" for further information',
                                  '''EGFR, Q9Y6I3 - 1, ENSG00000074410, ENSP00000075430.7, HGNC: 10728, UPI00022F85F1, NP_004702.2, SEMA5B - 211, BRAF_HUMAN, NM_001304833''',key=ss.run_id)
        else:
            input1 = st.text_area('Paste any Ensembl/Uniprot/Refseq ID\'s, gene names or a raw amino acid sequence: ', '''''',key=ss.run_id)
        file_upload, search_button = st.beta_columns([2.58,1])
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
            with st.spinner('Checking database . . .'):
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
            st.markdown("### Alignment Preview")
            st.write('\n')
            chosen_gene = list(nested_dict.keys())[0]
            index_gene_object = list(list(nested_dict.values())[0].keys())[0]
            transcript_list, index_gene = Visualise_Alignment.fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects, nested_dict, chosen_gene)
            st.write('Number of Isoform Entries for '+chosen_gene+':',len(transcript_list))
            reference_select, whitespace = st.beta_columns([1, 2.2])
            with reference_select:
                chosen_reference = st.selectbox('Choose your reference transcript: ',[transcript[0] for transcript in transcript_list])
                index_of_reference_transcript = Visualise_Alignment.get_index_of_chosen_transcript(chosen_reference,transcript_list)
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, index_gene_object)
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
            st.write('\n')
            st.text('\n')
            with st.spinner('Visualising Alignments . . .'):
                Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript,list_of_gene_objects,index_gene_object,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            #Table section
            parameter_change = False
            chosen_columns = Input_flow.chose_columns(nested_dict,dict_of_IDs,ss.run_id_table,parameter_change)
            generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,list_of_gene_objects,index_gene_object,chosen_columns,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length=exon_length_AA)
            slot1 = st.empty()
            value = Table_Generation.display_filter_option_AA()
            if value == "":
                    slot1.write(generated_table)
            else:
                filter_df = Table_Generation.filter_all_columns_of_df(value, generated_table)
                if not filter_df.empty:
                    slot1.write(filter_df)
                    st.info('ℹ️ Delete value to go back to original mapping table.')
                else:
                    st.warning('Value "'+str(value)+'" does not exist in the dataframe.')
            Input_flow.generate_download_section(generated_table)


        #case of using multiple ID's
        elif ss.searched_clicked and len(input1_IDs) > 1 and not no_elements:
            using_IDs = True
            st.markdown("### Alignment Preview")
            st.text('\n')
            genes, reference = st.beta_columns([2,1.7])
            with genes:
                chosen_gene = st.selectbox('Select Gene:',Visualise_Alignment.create_list_gene_selection(list_of_gene_objects,nested_dict))
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
            with st.spinner('Visualising Alignments . . .'):
                Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript, list_of_gene_objects,gene_index, match, mismatch, open_gap_penalty, gap_extension_penalty,exon_length_AA)
            # Table section
            parameter_change=False
            if  [match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA] != ss.parameters:
                parameter_change = True
                ss.parameters = [match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA]
            chosen_columns = Input_flow.chose_columns(nested_dict,dict_of_IDs,ss.run_id_table,parameter_change)
            if chosen_columns:
                df_all = Table_Generation.create_table_for_dict_of_gene_objects(nested_dict,list_of_gene_objects,chosen_columns, match, mismatch, open_gap_penalty, gap_extension_penalty)
                if not df_all.empty:
                    with st.spinner('Preparing Preview of Mapping Table . . .'):
                        slot1 = st.empty()
                        value = Table_Generation.display_filter_option_AA()
                        if value=="":
                            slot1.write(df_all)
                        else:
                            filter_df = Table_Generation.filter_all_columns_of_df(value, df_all)
                            if not filter_df.empty:
                                slot1.write(filter_df)
                                st.info('ℹ️ Delete value to go back to original mapping table.')
                            else:
                                st.warning('Value "' + str(
                                    value) + '" does not exist in the dataframe.')
                    st.text('\n')
                    Input_flow.generate_download_section(df_all)
                else:
                    st.warning('No amino acid positions mapped currently. Tweak function parameters to generate matches.')
            ss.run_id_table += 1

        #Input 2 Area
        if using_IDs== False and raw_aa:
            input2 = st.text_area('Paste raw amino acid sequence of alternative isoform and click align: ', '''''', key=ss.run_id)
            align=st.button('Align')
            if align:
                ss.align_clicked = True
                ss.searched_clicked = False
            if input1 != "" and input2 != "" and ss.align_clicked and ss.searched_clicked==False:
                #Sidebar pop up, make function out of it?
                match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, raw=True)
                st.markdown("### Results")
                #st.write("\n")
                #st.markdown("##### Unfiltered Alignment:")
                #st.write("\n")
                needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(input1,input2, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)
                isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
                #st.text(Alignment_preview)
                st.write("\n")
                st.markdown("#### Alignment")
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
        placehold, clear_all = st.beta_columns([6.1, 1])
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
        st.write('The Restful API is accessible trough the url www.isoaligner.org/api. Currently, a get method called /map for the retrieval of mapping tables between corresponding amino acid position position is available as well as the method /align to retrieve the alignment of two raw protein sequences.')

        st.write("--------------------------")
        st.header("Downloads")
        st.write('soon available')
        st.write("--------------------------")

    if choice == 'Manual & About':

        st.title(" Amino Acid Isoform Aligner")
        st.markdown("### Biologically Appropriate Alignment of Isoforms:")
        st.write('The challenge of aligning protein sequences of two isoforms is essentially matching the exons correctly.'
                 ' The IsoAligner algorithm exploits the biological characteristics of isoforms by fitting the parameters of the Needleman Wunsch global alignment and validating its solution in an extra step to assure positional mapping of biologically corresponding amino acid positions only.  ')
        problem_schema = Image.open('Pictures/mapping_problem_may.png')
        st.image(problem_schema,use_column_width=True)
        st.write('To avoid and discard falsely mapped positions of distinct exons (e.g. Exon4 and Exon5) the parameters of the alignment are tweaked as follows:')
        needleman, minimal = st.beta_columns([1, 1.])
        with needleman:
            st.markdown("#### Needleman-Wunsch global alignment")
            st.write(" - Big open gap penalty (Default -1.75) \n"
                 "- Small extend gap penalty (Default 0)\n"
                    "- Normal match and mismatch values (Default 1, -2)" )
        with minimal:
            st.markdown("#### Discard falsely matched positions")
            st.write('- By definition of a minimal exon length in numbers of consecutive AA. The length is gene-specific or at least 5 AA per default.')
        st.markdown("### Alignment example:")
        st.write("First off, IsoAligner aims at exon pattern alignment solutions. The generated AA matches are then additionally validated by the minimal exon length function. Alignment sections only containing partial diffuse mapping are being recognised as random matches and are marked as 'x' and ultimately discarded.")
        example = Image.open('Pictures/example_may.png')
        st.image(example, use_column_width=True)
        st.write("--------------------------")
        st.markdown("### Manual Alignment Tool")
        st.write("Quick Start: Click on 'Show Example' and then 'Search Library for IDs' to get a overview.")
        st.write("1. Paste ID's, gene names or raw amino acid sequences"
                 "\n    - The current human library consists of ~18K protein coding genes covering +130K protein sequences and +1.3M mapped Isoforms ID's from Ensembl, Uniprot, Refseq & HGNC. Included are Gene, Transcript & Protein ID's of various types (see figure below)."
                 "\n    - Click 'Search Library for IDs' or 'Align' to compute alignments"
                 "\n 2. Tweak function parameters in the sidebar and inspect the alignment previews"
                 "\n    - Set the mininmal exon length (in AA)"
                 "\n""    - Set Needleman-Wunsch parameters"
                 "\n 3. Explore the computed mapping table."
                "\n    - Filter for specific amino acid position "
                "\n     - Select ID's to be included dataframe."
                "\n    - Download dataframe: tab or comma separated (tsv/csv) "
                 )
        st.markdown("#### ⚠️ Important:")
        st.write("When multiple ID's are entered, the reference transcript for the generation of the mapping table is automatically chosen, unless a specific transcript or protein ID is used in the input field."
                 " In this case, the isoform with the longest sequence is used as the reference to align against."
                 " Also, be aware that the minimal exon length for the generation of the mapping table is likewise automatically chosen in this context. (gene-specific or at least 5 AA).")
        st.write("--------------------------")
        st.markdown("### Human Isoform Library Overview:")
        st.write('\n')
        total_number_of_genes, total_number_of_isoforms, genes_without_isoforms,Ids_in_total, minimal_exon_lengths, mean_exon = Statistics.list_of_gene_objects_statistics(list_of_gene_objects)
        picture, statistics = st.beta_columns([2, 1.2])
        with statistics:
            st.markdown("#### Statistics:")
            st.write('\n')
            st.write('Genes: ',total_number_of_genes)
            st.write('Isoforms in total: ',total_number_of_isoforms)
            st.write(" ID's in total:", Ids_in_total)
            st.write('Ø Isoform per gene:', round(total_number_of_isoforms/total_number_of_genes,1))
            st.write("Ø ID's per AA seq:", round(Ids_in_total/total_number_of_isoforms,1))
            st.write('Genes with minimal exons:', minimal_exon_lengths)
            st.write('Median minimal exon length:', mean_exon)
        with picture:
            st.markdown("#### Generation and Structure of Library:")
            st.write('\n')
            library = Image.open('Pictures/human_libary.png')
            st.image(library, use_column_width=True,width=None)
        #st.write("Gene object attributes:", Gene.list_of_attributes())
        #st.write("Collection of Isoform ID's includes:",Protein_isoform.list_of_attributes())
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
        st.write('---------------')
        st.markdown("#### License:")
        st.write("\n")
        html_string = '<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">IsoAligner</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://www.linkedin.com/in/jacob-hanimann-778032137/" property="cc:attributionName" rel="cc:attributionURL">Jacob Hanimann</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.'
        st.markdown(html_string, unsafe_allow_html=True)
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

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

#declare session state variables
ss = SessionState.get(clicked=False,searched_clicked=False, align_clicked=False, generate=False,run_id=0,example=False, clear_button=False)

#import database
list_of_gene_objects = Input_flow.import_data_from_github('https://drive.google.com/file/d/18Rs-GIQgU06AA2yJyB5-g4ctmqqGyEm8/view?usp=sharing')

#new link = https://drive.google.com/file/d/18Rs-GIQgU06AA2yJyB5-g4ctmqqGyEm8/view?usp=sharing
# old link = https://drive.google.com/file/d/1hKq8ksU1I8A61zvqOkKrx8WIfbP32TG4/view?usp=sharing

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
        header, tRNA = st.beta_columns([3, 1.3])
        with header:
            # Title
            st.title(" Amino Acid Isoform Aligner")
            st.subheader("Map amino acid positions across isoforms")
            st.write(
                "Align isoforms with the Needleman-Wunsch algorithm and set the minimal exon length to discard falsely mapped positions (random matches) of two distinct exons."
                " The table of correctly mapped positions can be downloaded as a file in several formats. A preview of the alignments is displayed dynamically. ")
        with tRNA:
            st.write('\n')
            st.write('\n')
            image2 = Image.open('TransferMessengerRNA.tif')
            st.image(image2,use_column_width=True, caption='TransferMessengerRNA & Ribosome in Action')


        st.write("--------------------------")
        st.sidebar.markdown("### 🧬️Organism")
        st.sidebar.selectbox('Select species', ['🧍🏽Homo Sapiens', '🐁 Mouse', '🪰 D. Melanogaster',])
        st.sidebar.write("--------------------------")

        #fixed in put area
        title, example_button = st.beta_columns([4.1,1])
        with title:
            st.markdown("#### Input")
        with example_button:
            if st.button('Show Example'):
                ss.example = True
                ss.searched_clicked = False

        if ss.example:
            input1 = st.text_area('Paste multiple ID\'s comma or newline and click on search library for ID\'s. Go to "Manual" for further information', '''EGFR, KRAS, Q9Y6I3, FUCA2, CD9, ENSG00000074410, FAM168A, ENSP00000075430.7''',key=ss.run_id)
        else:
            input1 = st.text_area('Paste any Ensembl/Uniprot/Refseq ID\'s, gene names or a raw amino acid sequence: ', '''''',key=ss.run_id)
        file_upload, search_button = st.beta_columns([2.4,1])
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

        #case of using one ID
        if ss.searched_clicked and bool(input1_IDs) and len(input1_IDs) == 1 and list(input1_IDs.values())[0] != 'not found': #check if dictionary is not empty
            using_IDs = True
            #st.write(input1_IDs)
            #st.write(list(input1_IDs.values())[0])
            st.markdown("### Alignments")
            reference_select, placeholder = st.beta_columns([1,2.5])
            with reference_select:
                chosen_gene = list(nested_dict.keys())[0]
                index_gene_object = list(list(nested_dict.values())[0].keys())[0]
                transcript_list = Visualise_Alignment.fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,nested_dict,chosen_gene)
                chosen_reference = st.selectbox('Choose your reference transcript: ',[transcript[0] for transcript in transcript_list])
                index_of_reference_transcript = Visualise_Alignment.get_index_of_chosen_transcript(chosen_reference,transcript_list)
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters()
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
            st.write('\n')
            st.text('\n')
            Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript,list_of_gene_objects,index_gene_object,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            #Table section
            chosen_columns = Input_flow.chose_columns()
            generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,list_of_gene_objects,index_gene_object,chosen_columns,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            st.text('\n')
            st.write(generated_table)
            st.text('\n')
            Input_flow.generate_download_section(generated_table)


        #case of using multiple ID's
        elif ss.searched_clicked and len(input1_IDs) > 1:
            using_IDs = True
            st.markdown("### Alignments")
            st.text('\n')
            genes, reference = st.beta_columns([2,1.7])
            with genes:
                chosen_gene = st.selectbox('Select Gene',Visualise_Alignment.create_list_gene_selection(list_of_gene_objects,nested_dict))
            with reference:
                transcript_list = Visualise_Alignment.fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,nested_dict, chosen_gene)
                chosen_reference = st.selectbox('Choose your reference transcript: ',[transcript[0] for transcript in transcript_list])
                index_of_reference_transcript = Visualise_Alignment.get_index_of_chosen_transcript(chosen_reference,transcript_list)
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters()
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
            df_all = Table_Generation.create_table_for_dict_of_gene_objects(nested_dict,list_of_gene_objects,chosen_columns, match, mismatch, open_gap_penalty, gap_extension_penalty,exon_length_AA)
            if not df_all.empty:
                st.write(df_all)
                st.text('\n')
                Input_flow.generate_download_section(df_all)
            else:
                st.warning('No amino acid positions mapped currently. Tweak function parameters to generate matches.')


        #Input 2 Area
        if using_IDs== False and raw_aa:
            input2 = st.text_area('Paste Amino Acid sequence of alternative isoform: ', '''''', key=ss.run_id)
            align=st.button('Align')
            if align:
                ss.align_clicked = True
                ss.searched_clicked = False
            if input1 != "" and input2 != "" and ss.align_clicked and ss.searched_clicked==False:
                #Sidebar pop up, make function out of it?
                match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters()
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
        st.write("--------------------------")
        st.header("Downloads")
        st.write("--------------------------")

    elif choice == 'Manual & About':
        st.title(" Amino Acid Isoform Aligner")
        st.write("--------------------------")
        st.markdown("#### Problem Statement:")
        image2 = Image.open('Problem_statement.png')
        st.image(image2,use_column_width=True)
        st.write("\n")
        st.markdown("#### About the code:")
        st.write("\n")
        st.markdown("##### Needleman-Wunsch Algorithm:")
        st.write("It is also sometimes referred to as the optimal matching algorithm and the global alignment technique. The Needleman–Wunsch algorithm is still widely used for optimal global alignment, particularly when the quality of the global alignment is of the utmost importance.")
        st.markdown("##### check_for_wrong_exon_alignments function :")
        st.write(" This function helps to identify falsely aligned elements (distinct exons) when globally aligning isoforms of a gene. The Needleman Wunsch algorithm randomly alignes fractions of non-identical exons since the optimization of the algorithm is to maximize matches (which makes sense with homologues but not with isoforms). This function discards such fractions of the alignment by rejecting exons shorter than the defined minimal length (in AA).")
        st.write("\n")
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
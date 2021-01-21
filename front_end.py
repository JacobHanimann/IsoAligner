import pandas as pd
import streamlit as st
from functions import * #import all functions from the functions.py file
import pickle
import sys
import SessionState

#declare session state variables
ss = SessionState.get(clicked=False,searched_clicked=False, align_clicked=False, generate=False,run_id=0)


#move classes from database to functions script






@st.cache(allow_output_mutation=True)
def import_data(file):
    with open(file,"rb") as fp:  # Pickling
           list_of_gene_objects = pickle.load(fp)
    return list_of_gene_objects


#fetch data from google drive
#url = 'https://drive.google.com/file/d/1C1To0a_y88LG811fIqvivvB4dT2SNF8G/view?usp=sharing'
#st.write(url.split('/')[-2])
#path = 'https://drive.google.com/uc?export=download&id='+url.split('/')[-2]
#st.write('https://drive.google.com/uc?export=download&id='+url.split('/')[-2])
#df = pd.read_csv(path)

list_of_gene_objects = import_data('list_of_gene_objects_with_fasta.txt')


#Playground
#if st.button("Search"):
#    ss.clicked = True
#if st.button("Reset"):
#    ss.clicked = False
#
#st.write(ss.clicked)
#
#options = st.multiselect(
#    'What are your favorite colors',
#    ['Green', 'Yellow', 'Red', 'Blue'],
#    ['Yellow', 'Red'])
#
#st.write('You selected:', options)
#e = RuntimeError('This is an exception of type RuntimeError')
#st.exception(e)
#st.success('This is a success message!')
#st.info("This is information")
#
#
#first_list=split_elements_from_user_input_string('ENSG00000282353,ENSG00000003137,ENSG00000006606,ENSG00000003137.8')
#st.write(first_list)
#
#second_dict=identify_IDs_from_user_text_input('ENSG00000282353,ENSG00000003137,ENSG00000006606,ENSG00000003137.8, KRAS, ENST00000004531,ENSP00000005178.5, Q96HP8,UPI00001BDC11,O14792,Q9Y258')
#st.write(second_dict)
#index= search_through_database_with_known_ID_Type(list_of_gene_objects,identify_IDs_from_user_text_input('ENSG00000282357\nENSG00000003137\nENSG00000006606\nENSG00000003137.8\nENSG00000003509.16\nHBB'))
#st.write(index)
#
#color = st.select_slider(
#'Select a color of the rainbow',
#options=['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet'])
#st.write('My favorite color is', color)
#
#
#title = st.text_input('Movie title', 'Life of Brian')
#st.write('The current movie title is', title)



#image = Image.open('sunrise.jpg')
#st.image(image, caption='Sunrise by the mountains',
         #use_column_width=True)

#Streamlit website

def main():
    """ Isoform Alignment Tool """
    #Title
    st.title(" Amino Acid Isoform Aligner")

    #Sidebar
    activity = ['Alignment Tool', 'Download Pre-Computed Data', 'About & Source Code']
    st.sidebar.markdown("## Navigation")
    choice = st.sidebar.radio("Go to", activity)
    st.sidebar.write("\n")


    #Alignment tool section
    if choice == 'Alignment Tool':
        #image1 = Image.open('Isoform_picture.png')
        #st.image(image1,use_column_width=True)
        st.subheader("A simple tool to map amino acid positions across isoforms")
        st.write("Align isoforms with the Needleman-Wunsch algorithm and set the minimal exon length to discard falsely mapped positions (random matches) of two distinct exons."
                 " The table of correctly mapped positions can be downloaded as a file in several formats. A preview of the alignments is displayed dynamically. ")

        st.write("--------------------------")
        st.sidebar.markdown("### 🧬️Organism")
        st.sidebar.selectbox('Select species', ['Homo Sapiens 🧍🏽‍', 'D. Melanogaster 🪰', 'Mouse 🐁', 'Frog 🐸', 'Mermaid 🧜🏼‍'])
        st.sidebar.write("--------------------------")

        #fixed in put area
        title, example_button = st.beta_columns([4.1,1])
        with title:
            st.markdown("#### Input")
        with example_button:
            st.button('Load Example')
        input1 = st.text_area('Paste gene names, IDs or raw amino acid sequence of reference isoform: ', '''EGFR''',key=ss.run_id)
        input1_IDs =  search_through_database_with_known_ID_Type(list_of_gene_objects,identify_IDs_from_user_text_input(input1))
        file_upload, search_button = st.beta_columns([2.4,1])
        with file_upload:
            agree = st.checkbox("Click here to upload list of gene names or ID's")
            if agree:
                input1 = st.file_uploader("Accepted ID's: Ensembl, Refseq, Uniprot (Accession/Uniparc)", type=["gz", "txt"])
        with search_button:
            search = st.button('Search Database for IDs')
            if search:
                ss.searched_clicked =True

        #set default for displaying second text_area input for input2
        using_IDs= False

        #check what user input is

        #case of using one ID's
        if ss.searched_clicked and bool(input1_IDs) and len(input1_IDs) == 1 and list(input1_IDs.values())[0] != 'not found': #check if dictionary is not empty
            using_IDs = True
            #st.write(input1_IDs)
            #st.write(list(input1_IDs.values())[0])
            st.markdown("### Alignments")
            reference_select, placeholder = st.beta_columns([1,2.5])
            with reference_select:
                chosen_reference = st.selectbox('Choose your reference transcript: ',fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,list(input1_IDs.values())[0]))
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = sidebar_pop_up_parameters()
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.write('\n')
            st.text('\n')
            display_alignment_for_one_gene_from_database(chosen_reference,list_of_gene_objects,list(input1_IDs.values())[0],match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            st.markdown("#### Table")
            options = st.multiselect(
                'Choose columns',
                ['Gene name', 'Transcript ID', 'AA', 'Ref Position', 'Isoform Position','Refseq ID', 'Uniprot_ID'],
                ['Gene name', 'Transcript ID', 'AA', 'Ref Position'])


        #case of using multiple ID's
        elif ss.searched_clicked and len(input1_IDs) > 1:
            st.write('multiple IDs')
            using_IDs = True
            st.write(input1_IDs)
            st.write(list(input1_IDs.values())[0])
        #case user types in aminoacid and clicks on search database
        elif ss.searched_clicked and extract_only_AA_of_Fasta_file(input1)!=None and ss.align_clicked==False:
            st.warning("Looks like an Amino Acid sequence! Paste in your second sequence below and click 'Align' ")
        elif ss.searched_clicked:
            st.warning("Couldn't find any ID's")
        st.write("\n")


        #Input 2 Area
        if using_IDs== False:
            input2 = st.text_area('Paste Amino Acid sequence of alternative isoform: ', '''''')
            align=st.button('Align')
            if align:
                ss.align_clicked = True
                ss.searched_clicked = False
            st.write("--------------------------")
            if input1 != "" and input2 != "" and ss.align_clicked and ss.searched_clicked==False:
                #Sidebar pop up, make function out of it?
                match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = sidebar_pop_up_parameters()
                st.markdown("### Results")
                #st.write("\n")
                #st.markdown("##### Unfiltered Alignment:")
                #st.write("\n")
                maped_tuple = map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(input1, input2, match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
                #st.text(Alignment_preview)
                st.write("\n")
                st.markdown("##### Alignment")
                st.write("\n")
                st.text(visualise_alignment_dynamically(maped_tuple[5],maped_tuple[6],maped_tuple[4]))
                st.markdown(" ###### ℹ️Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
                st.write("\n")
                st.write("\n")
                generated_table = create_pandas_dataframe_raw_aa_sequence(maped_tuple)
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
                    st.markdown(get_table_download_link(generated_table, 'dataframe.' + sep_choice, sep),unsafe_allow_html=True)
                with table:
                    st.markdown("##### Correctly mapped AA positions")
                    st.write("\n")
                    st.write(generated_table)
                st.write("--------------------------")

        #Clear all button
        placehold, clear_all = st.beta_columns([4.5, 1])
        with clear_all:
           reset= st.button('Clear All')
        if reset:
            ss.run_id +=1

    elif choice == 'Download Pre-Computed Data':
        st.header("Pre-Computed mapped isoforms")
        st.write("--------------------------")
        st.markdown("#### Refseq (4GB):")
        st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv','', 'Refseq_Isoforms.tsv'), unsafe_allow_html=True)
        st.markdown("#### Ensembl (4GB):")
        st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv','', 'Ensembl_Isoforms.tsv'), unsafe_allow_html=True)
        st.markdown("#### All ID's (8GB):")
        st.markdown(get_binary_file_downloader_html('/Users/jacob/Documents/GitHub/Mapping_Transcripts/streamlitmapping.tsv','', 'All_Isoforms.tsv'), unsafe_allow_html=True)

    elif choice == 'About & Source Code':
        st.write("--------------------------")
        st.markdown("#### Problem Statement:")
        image2 = Image.open('Problem_statement.png')
        st.image(image2,
        use_column_width=True)
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
        #st.text('LinkedIn: ')
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

#with open("list_of_gene_objects_with_fasta.txt", "rb") as fp:   #Pickling
    #list_of_gene_objects_with_fasta = pickle.load(fp)

#for gene in list_of_gene_objects_with_fasta:
    #if len(gene.protein_sequence_isoform_collection) >0:
        #print(gene.gene_symbol)


#Execution
if __name__ == '__main__':
    main()


#for gene in list_of_gene_objects_with_fasta:
#    if type(gene.protein_sequence_isoform_collection) == list:
#        if len(gene.protein_sequence_isoform_collection) >0:
#            print(gene.ensembl_gene_symbol)

#for gene in list_of_gene_objects:
#    print(gene.ensembl_gene_symbol)
#    print(gene.previous_symbols)
#    #print(type(gene.alias_symbols),type(gene.previous_symbols))
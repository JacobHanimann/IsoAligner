import pandas as pd
import streamlit as st
from functions_old import * #import all functions from the functions_old.py file
import pickle
import sys
import SessionState
from Streamlit_community import *

#declare session state variables
ss = SessionState.get(clicked=False,searched_clicked=False, align_clicked=False, generate=False,run_id=0,example=False, clear_button=False)


#move classes from database to functions script

#st.write('session state')
#st.write(ss.run_id)


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

list_of_gene_objects = import_data('../list_of_gene_objects_with_fasta.txt')


#Playground

#st.slider("Slide me!", 0, 100, key=ss.run_id)
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
        if ss.searched_clicked ==False:
            with example_button:
                if st.button('Show Example'):
                    ss.example = True
                    ss.searched_clicked = False

        if ss.example:
            input1 = st.text_area('Paste gene names, IDs or raw amino acid sequence of reference isoform: ', '''EGFR, KRAS, Q9Y6I3, FUCA2, CD9, ENSG00000074410, FAM168A, ENSP00000075430.7''',key=ss.run_id)
        else:
            input1 = st.text_area('Paste gene names, IDs or raw amino acid sequence of reference isoform: ', '''''',key=ss.run_id)
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

        if ss.searched_clicked:
            dict_of_IDs = identify_IDs_from_user_text_input(input1)
            input1_IDs = search_through_database_with_known_ID_Type(list_of_gene_objects, dict_of_IDs)
            nested_dict = generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs, input1_IDs,
                                                                                            list_of_gene_objects)
            #nested_dict = remove_dict_elements_with_no_gene_object_match(nested_dict)
        #check what user input is

        #case of using one ID's
        if ss.searched_clicked and bool(input1_IDs) and len(input1_IDs) == 1 and list(input1_IDs.values())[0] != 'not found': #check if dictionary is not empty
            using_IDs = True
            #st.write(input1_IDs)
            #st.write(list(input1_IDs.values())[0])
            st.markdown("### Alignments")
            reference_select, placeholder = st.beta_columns([1,2.5])
            with reference_select:
                chosen_gene = list(input1_IDs.keys())[0]
                index_gene_object = list(list(input1_IDs.values())[0].keys())[0]
                chosen_reference = st.selectbox('Choose your reference transcript: ',fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,nested_dict,chosen_gene))
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = sidebar_pop_up_parameters()
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
            st.write('\n')
            st.text('\n')
            display_alignment_for_one_gene_from_database(chosen_reference,list_of_gene_objects,index_gene_object,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            st.markdown("#### Mapped Amino Acid Positions Table")
            chosen_columns = st.multiselect(
                'Select columns',
                ['Gene name', 'Ensembl Gene ID','Ensembl Transcript ID','Ensembl Protein ID','Refseq Gene ID','Refseq Transcript ID','Uniprot Accession ID','Uniprot Isoform ID', 'Uniparc ID','Ensembl Gene ID version', 'Ensembl Transcript ID version', 'Ensembl Protein ID version','HGNC gene symbol'],
                ['Gene name', 'Ensembl Protein ID'])
            generated_table = create_table_for_one_gene_object(chosen_reference,list_of_gene_objects,index_gene_object,chosen_columns,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
            st.text('\n')
            st.write(generated_table)
            st.text('\n')
            download, format = st.beta_columns([0.19,1])
            with download:
                st.markdown("#### 📁 Download")
                sep_choice = st.radio('Choose file format:', ['tsv', 'csv'])
                if sep_choice == "tsv":
                    sep = '\t'
                else:
                    sep = ','
            with format:
                st.text('\n')
                st.text('\n')
                st.text('\n')
                st.text('\n')
                st.markdown(get_table_download_link(generated_table, 'DataframeMappedIsoforms.' + sep_choice, sep),unsafe_allow_html=True)


        #case of using multiple ID's
        elif ss.searched_clicked and len(input1_IDs) > 1:
            using_IDs = True
            st.markdown("### Alignments")
            st.text('\n')
            genes, reference = st.beta_columns([1,1.25])
            with genes:
                chosen_gene = st.selectbox('Select Gene',[element+' ('+str(len(list_of_gene_objects[list(index.keys())[0]].protein_sequence_isoform_collection))+' Isoforms)' for element,index in input1_IDs.items()])
            with reference:
                chosen_reference = st.selectbox('Select reference isoform', fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects,nested_dict,chosen_gene))
            ss.generate = True
            st.text('\n')
            match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = sidebar_pop_up_parameters()
            st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
            st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
            st.write('\n')
            st.text('\n')
            gene_index = list(nested_dict[re.split(' \(',chosen_gene)[0]])[0]
            #st.write('indexes of gene objects:')
            #st.write(nested_dict)
            display_alignment_for_one_gene_from_database(chosen_reference, list_of_gene_objects,gene_index, match, mismatch, open_gap_penalty, gap_extension_penalty,exon_length_AA)
            st.markdown("#### Mapped Amino Acid Positions Table")
            chosen_columns = st.multiselect('Select further columns',['Gene name', 'Ensembl Gene ID', 'Ensembl Transcript ID', 'Ensembl Protein ID', 'Refseq Gene ID', 'Refseq Transcript ID', 'Uniprot Accession ID', 'Uniprot Isoform ID', 'Uniparc ID', 'Ensembl Gene ID version', 'Ensembl Transcript ID version', 'Ensembl Protein ID version', 'HGNC gene symbol'],['Gene name', 'Ensembl Protein ID'])
            df_all = create_table_for_dict_of_gene_objects(nested_dict,list_of_gene_objects,chosen_columns, match, mismatch, open_gap_penalty, gap_extension_penalty,exon_length_AA)
            st.write(df_all)
            st.text('\n')
            download, format = st.beta_columns([0.19, 1])
            with download:
                st.markdown("#### 📁 Download")
                sep_choice = st.radio('Choose file format:', ['tsv', 'csv'])
                if sep_choice == "tsv":
                    sep = '\t'
                else:
                    sep = ','
            with format:
                st.text('\n')
                st.text('\n')
                st.text('\n')
                st.text('\n')
                st.markdown(get_table_download_link(df_all, 'DataframeMappedIsoforms.' + sep_choice, sep),
                            unsafe_allow_html=True)


        #case user types in aminoacid and clicks on search database
        elif ss.searched_clicked and extract_only_AA_of_Fasta_file(input1)!=None and ss.align_clicked==False:
            st.warning("Looks like an Amino Acid sequence! Paste in your second sequence below and click 'Align' ")
        elif ss.searched_clicked:
            st.warning("Couldn't find any ID's")
        st.write("\n")


        #Input 2 Area
        if using_IDs== False:
            input2 = st.text_area('Paste Amino Acid sequence of alternative_ID isoform: ', '''''', key=ss.run_id)
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
                needleman_mapped = map_FMI_on_COSMIC_Needleman_Wunsch_with_exon_check(input1,input2, match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA)
                isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
                #st.text(Alignment_preview)
                st.write("\n")
                st.markdown("##### Alignment")
                st.write("\n")
                percentage_reference, percentage_isoform = calculate_percentage_of_mapped_positions(isoform_pattern_check,input1,input2)
                st.text(visualise_alignment_dynamically(alignment_reference_fasta,alignment_isoform_fasta,isoform_pattern_check,percentage_reference,percentage_isoform))
                st.write("\n")
                st.markdown(" ###### ℹ️Syntax: 'x' are discarded matches determined by the minimal exon length and '|' are valid matches of identical exons")
                st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per isoform")
                st.write("\n")
                st.write("\n")
                generated_table = create_pandas_dataframe_raw_aa_sequence(needleman_mapped)
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

        #Clear all button
        st.write("--------------------------")
        placehold, clear_all = st.beta_columns([4.5, 1])
        with clear_all:
           if st.button('Clear All'):
              ss.clear_button = True
              ss.run_id +=1
              ss.example= False
              st.write('click')
              st.write(ss.run_id)
              ss.searched_clicked= False
              #raise st.script_runner.RerunException()



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
        image2 = Image.open('../Pictures/Problem_statement.png')
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
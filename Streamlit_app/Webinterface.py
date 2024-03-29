#Author: Jacob Hanimann
#Python 3.8

import sys
sys.path.append(".")
import SessionState
from Streamlit_Pop_ups import *
from IsoAligner_core.Visualise_Alignment import *
from IsoAligner_core.User_Input_Preparation import *
from IsoAligner_core.Input_flow import *
from IsoAligner_core.Table_Generation import *
from PIL import Image
from Statistics import *

#import database
list_of_gene_objects = Input_flow.import_data_from_github('Human_Isoform_Library/list_of_gene_objects_25th_july.txt.gz')

#declare session state variables
ss = SessionState.get(clicked=False, searched_clicked=False, align_clicked=False, generate=False, run_id=0, example=False, clear_button=False, run_id_table=1, parameters=[1, 2, 3, 4, 5], random_input = Input_flow.generate_random_example(list_of_gene_objects))

#Streamlit website
def main():
    """ Isoform Alignment Tool """

    #remove streamlit marks
    st.markdown(""" <style>
        footer {visibility: hidden;}
        </style> """, unsafe_allow_html=True)
    st.markdown("""
           <style>div[data-testid="stToolbar"] { display: none;}</style>
           """, unsafe_allow_html=True)

    #Sidebar
    activity = ['Alignment Tool', 'REST API & Downloads', 'Manual & About']
    st.sidebar.markdown("## Navigation")
    choice = st.sidebar.radio("Go to", activity)


    #Alignment tool section
    if choice == 'Alignment Tool':
        logo, name = st.columns([0.18,1])
        with logo:
            image2 = Image.open('Streamlit_app/Pictures/only_logo.png')
            st.image(image2, use_column_width=True, width=350)
        with name:
            st.markdown('# IsoAligner')
        st.header("Map Amino Acid Positions Across Isoforms")
        st.write(
            "Align protein isoforms interactively with a customized Needleman-Wunsch algorithm and gene-specific minimal exon lengths to retrieve clean positional mapping tables for corresponding amino acids in alternative splice variants."
            " The current human isoform library consists of ~19'000 protein coding genes covering ~120'000 protein sequences and +1.3 million mapped isoform IDs from Ensembl, UniProt, RefSeq and UCSC.")

        st.sidebar.markdown("### 🧬️Organism")
        st.sidebar.selectbox('Select species', ['🧍🏽Homo Sapiens'])

        #fixed in put area
        title, example_button = st.columns([2.3,1])
        with title:
            st.markdown("### Input")
            with example_button:
                if st.button('Show Random Example'):
                    ss.example = True

        if ss.example:
            input1 = st.text_area('Paste IDs (comma or newline separated) and click on Search and Align. Go to "Manual" for further information.',ss.random_input,key=str(ss.run_id))
        else:
            input1 = st.text_area('Paste one, two (pairwise) or multiple Ensembl/UniProt/RefSeq IDs, gene names or a raw amino acid sequence: ', '''''',key=str(ss.run_id))
        file_upload, search_button = st.columns([3.35,1])
        with file_upload:
            file_wanted = st.checkbox("upload list of IDs or gene names")
            if file_wanted:
                input1 = st.file_uploader("Accepted IDs: Ensembl, Refseq, Uniprot, UCSC, HGNC (Accession/Isoform/Uniparc)", type=[ "txt"])
                if input1 is not None:
                    input1 = input1.getvalue().decode("UTF-8")
                    ss.searched_clicked = True
            raw_aa = st.checkbox("insert 2nd raw amino acid sequence manually")
        with search_button:
            search = st.button('Search and Align')
            if search:
                ss.searched_clicked =True

        #set default for displaying second text_area input for input2
        using_IDs= False
        no_elements = False
        input1_IDs = []
        #state variable
        mode_of_action = 'unknown'


        if ss.searched_clicked and input1 !="""""" and input1!=None:
            with st.spinner('Checking library . . .'):
                dict_of_IDs = Input_preparation.identify_IDs_from_user_text_input(input1)
                #st.write(dict_of_IDs)
                input1_IDs = Input_flow.search_through_database_with_known_ID_Type(list_of_gene_objects, dict_of_IDs)
                #st.write(input1_IDs)
                Input_flow.show_which_elements_were_not_found(input1_IDs)
                cleaned_input1_IDs=Input_flow.remove_dict_elements_with_no_gene_object_match(input1_IDs)
                nested_dict = Input_flow.generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs, cleaned_input1_IDs,list_of_gene_objects)
                # identify input as pairwise or one ID per gene or mixed and report
                mode_of_action = Input_flow.report_mode_of_action(nested_dict)
                if mode_of_action=='pairwise':
                    Input_flow.show_identical_elements(nested_dict, list_of_gene_objects)
                if len (nested_dict)==0 or mode_of_action=='stop':
                    no_elements = True


        #case of using one ID per gene
        if mode_of_action == "one_ID_per_gene":
            #case of using one ID overall
            if ss.searched_clicked and not no_elements and  len(input1_IDs) == 1:
                using_IDs = True
                st.write("---")
                title, clear_button = st.columns([5.9,1.])
                with title:
                    st.markdown("### Search Results")
                with clear_button:
                    st.write('\n')
                    if st.button('Clear All'):
                        ss.clear_button = True
                        ss.random_input = Input_flow.generate_random_example(list_of_gene_objects)
                        ss.run_id += 1
                        ss.example = False
                        ss.searched_clicked = False
                        Streamlit_community.rerun_script_from_top()
                chosen_gene = list(nested_dict.keys())[0]
                index_gene_object = list(list(nested_dict.values())[0].keys())[0]
                transcript_list, index_gene = Visualise_Alignment.fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects, nested_dict, chosen_gene, dict_of_IDs)
                one_isoform = False
                if len(transcript_list)==1:
                    one_isoform =True
                    Input_flow.pop_up_isoform_info(list_of_gene_objects, index_gene_object, one_isoform)

                if not one_isoform:
                    st.write('Number of Isoform Entries for '+chosen_gene+':',len(transcript_list))
                    reference_select, whitespace = st.columns([1, 0.9])
                    with reference_select:
                        chosen_reference = st.selectbox('Choose your reference isoform: ',[transcript[0] for transcript in transcript_list])
                        index_of_reference_transcript = Visualise_Alignment.get_index_of_chosen_transcript(chosen_reference,transcript_list)
                    with st.expander("View details about this Isoform Entry"):
                        Input_flow.pop_up_isoform_info(list_of_gene_objects, index_gene_object, one_isoform, index_of_reference_transcript)
                    ss.generate = True
                    st.markdown("##### Alignment Preview")
                    st.write("\n")
                    match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, index_gene_object)
                    with st.expander("View Alignment Visualisations"):
                        st.markdown("###### ℹ Syntax: 'x' are discarded matches and '|' are valid correspondences determined by the minimal exon length function. Uppercase 'X' display AA mismatches (comes with a warning message if occurring). The percentage score represents the ratio of correctly mapped positions over the total number of positions per sequence. The positions are annotated every 10th amino acid. AA position annotated in a gap region corresponds to the last prior AA in the sequence.")
                        with st.spinner('Visualising Alignments . . .'):
                            Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript,list_of_gene_objects,index_gene_object,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA)
                    #Table section
                    parameter_change = False
                    chosen_columns = Input_flow.chose_columns(list_of_gene_objects,nested_dict,dict_of_IDs,ss.run_id_table,parameter_change)
                    generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,list_of_gene_objects,index_gene_object,chosen_columns,match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length=exon_length_AA)
                    if not type(generated_table)==tuple:
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
                    else:
                        st.warning(
                            'No amino acid positions mapped currently.')
                        st.info(' Tweak function parameters on the left sidebar to allow matches!')

            #case of using one ID per gene
            #case of using multiple ID'
            elif ss.searched_clicked and len(input1_IDs) > 1 and not no_elements:
                using_IDs = True
                st.write("---")
                title, clear_button = st.columns([5.9, 1])
                with title:
                    st.markdown("### Search Results")
                with clear_button:
                    st.write('\n')
                    if st.button('Clear All'):
                        ss.clear_button = True
                        ss.run_id += 1
                        ss.example = False
                        ss.searched_clicked = False
                        Streamlit_community.rerun_script_from_top()
                Table_Generation.generate_multiple_IDs(nested_dict, list_of_gene_objects, dict_of_IDs, ss)

        #pairwise alignments
        elif mode_of_action == 'pairwise':
            dict_of_pairwise = Input_flow.create_dict_for_pairwise_mode(nested_dict, list_of_gene_objects)
            st.write("----")
            title, clear_button = st.columns([5.9, 1])
            with title:
                st.markdown("### Search Results")
            with clear_button:
                st.write('\n')
                if st.button('Clear All'):
                    ss.clear_button = True
                    ss.run_id += 1
                    ss.example = False
                    ss.searched_clicked = False
                    Streamlit_community.rerun_script_from_top()
            show_all_isoforms = st.checkbox('Show all associated Isoform Entries')
            if show_all_isoforms:
                for gene, elements in dict_of_pairwise.items():
                    nested_dict.pop(elements[1])
                Table_Generation.generate_multiple_IDs(nested_dict, list_of_gene_objects, dict_of_IDs, ss)
            else:
                genes, reference = st.columns([2, 2])
                with genes:
                    chosen_gene = st.selectbox('Select Gene:',
                                               Visualise_Alignment.create_list_gene_selection(list_of_gene_objects,
                                                                                              nested_dict, pairwise=True))
                with reference:
                    transcript_list = dict_of_pairwise[re.split(' \(', chosen_gene)[0]]
                    chosen_reference = st.selectbox('Choose your reference isoform: ', transcript_list)
                    index_of_reference_transcript = list(nested_dict[chosen_reference].values())[0]
                    for element in transcript_list:
                        if element == chosen_reference:
                            continue
                        else:
                            index_of_alternative_transcript = list(nested_dict[element].values())[0]
                    gene_index = list(nested_dict[chosen_reference].keys())[0]
                if len(transcript_list) == 1:
                    one_isoform = True
                    Input_flow.pop_up_isoform_info(list_of_gene_objects, one_isoform, gene_index)
                else:
                    one_isoform = False
                ss.generate = True
                with st.expander('Details about Gene and Isoform Entry'):
                    Input_flow.pop_up_isoform_info(list_of_gene_objects, gene_index, one_isoform,
                                                   index_of_reference_transcript)
                match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(
                    list_of_gene_objects, gene_index)
                st.markdown("##### Alignment Preview")
                st.write("\n")
                with st.expander("View Alignment Visualisations"):
                    st.markdown(
                        "###### ℹ Syntax: 'x' are discarded matches and '|' are valid correspondences determined by the minimal exon length function. Uppercase 'X' display AA mismatches (comes with a warning message if occurring). The percentage score represents the ratio of correctly mapped positions over the total number of positions per sequence. The positions are annotated every 10th amino acid. AA position annotated in a gap region corresponds to the last prior AA in the sequence.")
                    with st.spinner('Visualising Alignments . . .'):
                        transcript_list_w_index = [(element, list(nested_dict[element].values())[0]) for element in
                                                   transcript_list]
                        Visualise_Alignment.display_alignment_for_one_gene_from_database(index_of_reference_transcript,
                                                                                     list_of_gene_objects, gene_index,
                                                                                     match, mismatch, open_gap_penalty,
                                                                                     gap_extension_penalty, exon_length_AA,
                                                                                     pairwise=transcript_list_w_index,
                                                                                     streamlit=True)
                # Table section
                parameter_change = False
                if len(dict_of_pairwise) == 1:
                    chosen_columns = Input_flow.chose_columns(list_of_gene_objects, nested_dict, dict_of_IDs,
                                                              ss.run_id_table,
                                                              parameter_change)
                    df_all = Table_Generation.create_mapping_table_of_two_IDs(list_of_gene_objects, gene_index,
                                                                              index_of_reference_transcript,
                                                                              index_of_alternative_transcript,
                                                                              chosen_columns, match, mismatch,
                                                                              open_gap_penalty, gap_extension_penalty,
                                                                              exon_length_AA, two_ids=True)
                # multiple pairwise
                else:
                    if [match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA] != ss.parameters:
                        parameter_change = True
                        ss.parameters = [match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA]
                    chosen_columns = Input_flow.chose_columns(list_of_gene_objects, nested_dict, dict_of_IDs,
                                                              ss.run_id_table,
                                                              parameter_change)
                    if chosen_columns:
                        df_all = Table_Generation.create_mapping_table_of_two_IDs_dict(nested_dict, list_of_gene_objects,
                                                                                       chosen_columns, match, mismatch,
                                                                                       open_gap_penalty,
                                                                                       gap_extension_penalty)
                if not type(df_all) == tuple:
                    with st.spinner('Preparing Preview of Mapping Table . . .'):
                        slot1 = st.empty()
                        value = Table_Generation.display_filter_option_AA()
                        if value == "":
                            slot1.write(df_all)
                            # st.dataframe(df_all.style.highlight_(axis=0))
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
                    st.warning(
                        'No amino acid positions mapped currently.')
                    st.info(' Tweak function parameters on the left sidebar to allow matches!')
                ss.run_id_table += 1


        #Input 2 Area
        if using_IDs== False and raw_aa:
            input2 = st.text_area('Paste raw amino acid sequence of alternative isoform and click align: ', '''''', key=str(ss.run_id+1))
            align=st.button('Align')
            if align:
                ss.align_clicked = True
                ss.searched_clicked = False
            if input1 != "" and input2 != "" and ss.align_clicked and ss.searched_clicked==False:
                if len(input2)<7:
                    st.error('Protein sequence has to be at least 7 AA long. Please reload page.')
                match, mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA = Streamlit_pop_ups.sidebar_pop_up_parameters(list_of_gene_objects, raw=True)
                needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(input1.upper(),input2.upper(), match, mismatch, open_gap_penalty,gap_extension_penalty, exon_length_AA, streamlit=True)
                isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
                st.write("\n")
                title, clear_button = st.columns([5.5,1])
                with title:
                    st.markdown("### Alignment Preview")
                with clear_button:
                    st.write('\n')
                    if st.button('Clear All'):
                        ss.clear_button = True
                        ss.random_input = Input_flow.generate_random_example(list_of_gene_objects)
                        ss.run_id += 2
                        ss.example = False
                        ss.searched_clicked = False
                        Streamlit_community.rerun_script_from_top()
                st.write("\n")
                percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check,input1,input2)
                st.text(Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta,alignment_isoform_fasta,isoform_pattern_check,percentage_reference,percentage_isoform))
                st.write("\n")
                st.markdown(" ######  ℹ️ Syntax: 'x' are discarded matches and '|' are valid correspondences determined by the minimal exon length function")
                st.markdown(" ###### The percentage score represents the ratio of correctly mapped positions over the total number of positions per sequence")
                st.write("\n")
                st.write("\n")
                generated_table = Table_Generation.create_pandas_dataframe_raw_aa_sequence(needleman_mapped)
                if not generated_table.empty:
                    table, whitespace, download = st.columns([1.25,0.2,1])
                    with download:
                        st.write('\n')
                        st.write('\n')
                        st.write('\n')
                        Input_flow.generate_download_section(generated_table, no_column=True)
                    with table:
                        st.markdown("### Mapped Amino Acid Positions Table")
                        st.write("\n")
                        st.write(generated_table)
                else:
                    st.warning('No amino acid positions mapped currently. Tweak function parameters to generate matches.')


    elif choice == 'REST API & Downloads':
        logo, name = st.columns([0.18, 1])
        with logo:
            image2 = Image.open('Streamlit_app/Pictures/only_logo.png')
            st.image(image2, use_column_width=True, width=350)
        with name:
            st.markdown('# IsoAligner')
        st.header("REST API")
        st.write('The Restful API is accessible via the url www.isoaligner.org/api. Currently, a GET method called /map for the retrieval of mapping tables for IDs of the human isoform library is available, as well as the method /align to retrieve the alignment of two raw protein sequences.')
        st.write("#### Browse Function Features")
        resource,method, parameters = st.columns([1,1,1])
        with resource:
            resource = st.selectbox('Resource:', ['...org/api/map','...org/api/align'])
        with method:
            method = st.selectbox('Methods:', ['GET'])
        with parameters:
            if resource=="...org/api/map":
                parameter = st.selectbox('Parameters: ',['id1','id2','pos',"min_ex_len",'df_ids',"match","mismatch", "open_gap", "gap_ext", 'ℹ️ Show All Parameters'])
            else:
                parameter = st.selectbox('Parameters: ',['seq1','seq2',"min_ex_len", "match","mismatch", "open_gap", "gap_ext", 'ℹ️ show all parameters'])
        #default parameters
        match = 1
        mismatch = -2
        gap_extend = 0
        gap_open = -1
        exon_length_AA = 12
        if resource == "...org/api/map":
            st.markdown("### Resource: /map")
            st.write('With this resource, you can access the human isoform library, compute alignments with specified parametes and retrieve whole mapping tables in json format or just single AA positions. The only required parameter that must be sent in the request body is id1, respectively an Isoform ID of any kind.')
        if resource == "...org/api/align":
            st.markdown("### Resource: /align")
            st.write('With this resource, you can align two raw amino acid sequences sent with the request and retrieve a mapping table in json format. The required parameters are seq1 and seq2.')

        if parameter == 'id1' or parameter =='ℹ️ Show All Parameters':
            st.write(" ##### Parameter: id1")
            st.write('ID of any type (Ensembl, Refseq, Uniprot, UCSC) to access the isoforms of a gene of the human isoform library. '
                     'To define the reference protein sequence which all other splice variants are align against, use a specified splice variant identifier. Otherwise, the longest isoform is automatically chosen as the reference sequence. A list of the supported ID types can be found under Manual & About. As a simple request example:')
            st.markdown("request: <em> www.isoaligner.org/api/map?id1=EGFR-201</em> ",
                    unsafe_allow_html=True)

        if parameter == 'id2' or parameter =='ℹ️ Show All Parameters':
            st.write(" ##### Parameter: id2")
            st.write(
                'Specific isoform ID (Ensembl, Refseq, Uniprot, UCSC) to use as the alternative splice variant to align with the reference sequence of id1. A list of the supported ID types can be found under Manual & About. Example:')
            st.markdown("request: <em> www.isoaligner.org/api/map?id1=EGFR-201&id2=EGFR-207 </em> ",
                    unsafe_allow_html=True)
        if parameter == 'pos' or parameter == 'ℹ️ Show All Parameters':
            st.write(" ##### Parameter: pos")
            st.write("In case of stating two IDs in the request, you can retrieve single corresponding AA positions of the alternative isoform sequence. As an example:")
            st.markdown("request: <em> www.isoaligner.org/api/map?id1=EGFR-201&id2=EGFR-207&pos=1038 </em> ", unsafe_allow_html=True)
            st.markdown("response: <em> 993 </em> ", unsafe_allow_html=True)

        if resource == "...org/api/map":
            if parameter == 'min_ex_len' or parameter == 'ℹ️ Show All Parameters':
                st.write(" ##### Parameter: min_ex_len")
                st.write("The alignment parameter for the minimal exon length (consecutive AA) is gene-specific per default and can be manually defined followingly: ")
                st.markdown("request: <em> www.isoaligner.org/api/map?id1=EGFR-201&id2=EGFR-207&min_ex_len=23 </em> ",
                        unsafe_allow_html=True)

        if parameter == 'df_ids' or parameter == 'ℹ️ Show All Parameters':
            st.write(" ##### Parameter: df_ids")
            st.write(
                "Choose which sequence database IDs you want to include in the mapping table. Per default, the mapping table consists of the same type of IDs sent with the request. Available options are: [ensembl, refseq, uniprot, ucsc, hgnc]. Example to retrieve only ensembl and uniprot IDs:")
            st.markdown(
                "request: <em> www.isoaligner.org/api/map?id1=EGFR-201&id2=EGFR-207&df_ids=[ensembl,uniprot] </em> ",
                unsafe_allow_html=True)

        #align resource
        if parameter == 'seq1' or parameter =='seq2' or parameter == 'ℹ️ show all parameters':
            st.write(" ##### Parameters: seq1 and seq2")
            st.write("Reference and alternative raw amino acid sequences. Must be at least 7 AA's long, for example:")
            st.text('www.isoaligner.org/api/align?seq1=CRSSWTAAMELSAEYLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRHSLDPSMTIHCDMVITYGLDQLENCQTCGTDYIISVLNLLTLI&seq2=YLREKLQRDLEAEHVEVEDTTLNRCSCSFRVLVVSAKFEGKPLLQRH')

        if resource=="...org/api/align":
            if parameter == 'min_ex_len' or parameter == 'ℹ️ show all parameters':
                st.write(" ##### Parameter: min_ex_len")
                st.write(
                    "The alignment parameter for the minimal exon length (consecutive AA) is per default 12 and can be manually defined followingly: ")
                st.text('www.isoaligner.org/api/map?id1=EGFR-201&id2=EGFR-207&min_ex_len=23')

        if parameter == 'match' or parameter == 'ℹ️ Show All Parameters' or parameter == 'ℹ️ show all parameters':
            st.write(" ##### Parameter: match")
            st.write(
                "[Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) alignment parameter to reward matches. This value must be ≥ 0.")

        if parameter == 'mismatch' or parameter == 'ℹ️ Show All Parameters' or parameter == 'ℹ️ show all parameters':
            st.write(" ##### Parameter: mismatch")
            st.write(
                "[Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) alignment parameter to penalize mismatches. This value must be ≤ 0.")

        if parameter == 'open_gap' or parameter == 'ℹ️ Show All Parameters' or parameter == 'ℹ️ show all parameters':
            st.write(" ##### Parameter: open_gap")
            st.write(
                "[Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) alignment parameter to penalize opening a gap. This value must be ≤ 0.")

        if parameter == 'gap_extend' or parameter == 'ℹ️ Show All Parameters' or parameter == 'ℹ️ show all parameters':
            st.write(" ##### Parameter: gap_open")
            st.write(
                "[Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) alignment parameter to penalize extending a gap. This value must be ≤ 0.")

        st.write("--------------------------")
        st.header("Pre-Computed Mapped Human Isoform Library")
        st.write('Here you can download a dataframe containing the mapping tables of every gene that is part of the human isoform library. The isoform with the longest sequence is used as the reference to align against all other splice variants.')
        st.write('\n')
        total_number_of_genes, number_of_good_genes, total_number_of_isoforms, genes_without_isoforms,Ids_in_total, minimal_exon_lengths, mean_exon, two_isoform_number, Ids_two = Statistics.list_of_gene_objects_statistics(list_of_gene_objects)
        info,needleman, minimal = st.columns([1, 1,1])
        with info:
            st.markdown("##### Volume")
            st.write(' - Genes: ', number_of_good_genes," \n"
            ' - Isoform entries in total: ', two_isoform_number," \n"
            " - IDs in total:", Ids_two," \n"
            ' - Ø Isoform per gene:', round(two_isoform_number / number_of_good_genes, 1)," \n"
            " - Ø IDs per AA seq:", round(Ids_two / two_isoform_number, 1)," \n")

        with needleman:
            st.markdown("##### Needleman-Wunsch function")
            st.write(" - match: ",match," \n"
                     " - mismatch: ",mismatch," \n"
                     " - open gap penalty: ",gap_open," \n"
                     " - gap extend penalty: ", gap_extend, " \n")
        with minimal:
            st.markdown("##### Minimal exon length function")
            st.write('Gene-specific length for all genes. The median of the sum of all shortest exons is',mean_exon)

        st.write('### 📁 Download')
        st.write('The dataframe is ~385MB, tab-separated and zipped (.tsv.gz).')

        Streamlit_community.create_download_section_from_ext_link('1-djWhoz-Yadi1vVKI6AArfA65u8p2clz','Click here to start download')

    elif choice == 'Manual & About':
        logo, name = st.columns([0.18, 1])
        with logo:
            image2 = Image.open('Streamlit_app/Pictures/only_logo.png')
            st.image(image2, use_column_width=True, width=350)
        with name:
            st.markdown('# IsoAligner')
        st.markdown("### Matching identical exons")
        st.write('The challenge of aligning protein sequences of two isoforms is essentially matching the exons correctly.'
                 ' The IsoAligner algorithm exploits the biological characteristics of isoforms by using the Needleman-Wunsch algorithm with an open gap penalty and validating its solution in an extra step to provide positional mapping of comprehensibly corresponding amino acids exclusively.  ')
        problem_schema = Image.open('Streamlit_app/Pictures/exon_problem.png')
        st.image(problem_schema,use_column_width=True)
        st.write('To avoid and discard falsely mapped positions of distinct exons (e.g. Exon4 and Exon5) the parameters of the alignment are tweaked as follows:')
        needleman, minimal = st.columns([1, 1.])
        with needleman:
            st.markdown("#### Needleman-Wunsch global alignment")
            st.write(" - Big open gap penalty (Default -1) \n"
                 "- No extend gap penalty (Default 0)\n"
                    "- Conventional match and mismatch values (Default 1, -2)" )
        with minimal:
            st.markdown("#### Discard falsely matched positions")
            st.write('- By definition of a minimal exon length in numbers of consecutive AAs. The length is gene-specific and at least 3.')
        st.markdown("### Alignment example:")
        st.write("First off, IsoAligner aims at exon pattern alignment solutions. The generated AA matches are then additionally validated by the minimal exon length function. Alignment sections only containing partial diffuse mapping are being recognised as random matches and are marked as 'x' and ultimately discarded from the mapping table.")
        example = Image.open('Streamlit_app/Pictures/example_june.png')
        st.image(example, use_column_width=True)
        st.write("--------------------------")
        st.markdown("### Manual Alignment Tool")
        st.write("Quick Start: Click on 'Show Example' and then 'Search and Align' to get a overview.")
        st.write("1. Enter either one isoform ID per gene **or** two isoform IDs per gene **or** a list of genes names **or** two raw amino acid sequences. The input can be tab, comma or whitespace separated."
                "\n   - Included are Gene, Transcript & Protein IDs of various types (see figure below). The current human library consists of ~19K protein coding genes covering ~120K protein sequences and ~1.3M mapped Isoforms IDs from Ensembl, Uniprot, Refseq & HGNC."                 "\n    - Click 'Search and Align' or 'Align' to compute alignments."
                 "\n 2. Tweak function parameters in the left sidebar and inspect the alignment previews."
                 "\n    - Set the mininmal exon length (in AA)."
                 "\n""    - Set Needleman-Wunsch parameters."
                 "\n 3. Explore the computed mapping table."
                "\n    - Filter for specific amino acid positions or transcript IDs."
                "\n     - Select ID types (Ensembl, UniProt, RefSeq, UCSC) to be included in the dataframe."
                "\n    - Download mapping table: tab or comma separated (tsv/csv). "
                 )
        st.markdown("#### ⚠️ Important:")
        st.write("When multiple IDs are entered from different genes, the reference transcript for the generation of the mapping table is automatically chosen, unless a specific transcript or protein ID is used in the input field."
                 " In this case, the isoform with the longest sequence is used as the reference to align against."
                 " Also, be aware that the minimal exon length for the generation of the mapping table is likewise automatically chosen in this context, gene-specific and at least 3 AA.")
        st.write("--------------------------")
        st.markdown("### Human Isoform Library Overview:")
        st.write('\n')
        total_number_of_genes, number_of_good_genes, total_number_of_isoforms, genes_without_mme,Ids_in_total, minimal_exon_lengths, mean_exon, two_isoform_number, Ids_two = Statistics.list_of_gene_objects_statistics(list_of_gene_objects)
        picture, statistics = st.columns([2, 1.2])
        with statistics:
            st.markdown("#### Statistics:")
            st.write('\n')
            st.write('Genes total:', minimal_exon_lengths)
            st.write('Isoform entries in total: ', total_number_of_isoforms)
            st.write('Genes with 2 ≤ entries: ',number_of_good_genes)
            st.write(" IDs in total:", Ids_in_total)
            st.write('Ø Isoform per gene:', round(total_number_of_isoforms/total_number_of_genes,1))
            st.write("Ø IDs per AA seq:", round(Ids_in_total/total_number_of_isoforms,1))
            st.write('Median minimal exon length:', mean_exon)
        with picture:
            st.markdown("#### Generation and Structure of Library:")
            st.write('\n')
            library = Image.open('Streamlit_app/Pictures/human_isoform_library.png')
            st.image(library, use_column_width=True,width=None)
        st.write("--------------------------")
        st.markdown("#### Contact:")
        st.write("\n")
        st.write("Please get in touch for suggestions or to report bugs :)")
        st.text('''Linkedin: Jacob Hanimann\nE-mail: jacobh@student.ethz.ch''')
        st.write('Bioinformatics group: https://www.sib.swiss/abdullah-kahraman-group')
        #Footnote
        st.write('---------------')
        st.markdown("#### License:")
        st.write("\n")
        html_string = '<a rel="license" href="http://creativecommons.org/licenses/by/4.0/" target="_blank" ><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" target="_blank" property="dct:title">IsoAligner</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://www.linkedin.com/in/jacob-hanimann-778032137/" target="_blank" property="cc:attributionName" rel="cc:attributionURL">Jacob Hanimann</a> is licensed under a <a rel="license" href="https://creativecommons.org/licenses/by/4.0/" target="_blank" >Creative Commons Attribution 4.0 International License </a>.'
        st.markdown(html_string, unsafe_allow_html=True)

#Execution
if __name__ == '__main__':
    main()

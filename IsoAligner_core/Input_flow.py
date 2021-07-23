import pickle
import urllib
from Streamlit_app.Streamlit_community import *
import gzip
import random
from IsoAligner_core.Gene import *
from IsoAligner_core.Protein_isoform import *

class Input_flow:
    pass

    @staticmethod
    @st.cache(allow_output_mutation=True)
    def import_data_from_url(url):
        path = 'https://drive.google.com/uc?export=download&id=' + url.split('/')[-2]
        file = urllib.request.urlopen(path)
        data = file.read()
        list_of_gene_objects = pickle.loads(data)
        return list_of_gene_objects


    @staticmethod
    #@st.cache(allow_output_mutation=True)
    def import_data_from_github(file):
        st.write(Gene())
        st.write('in the function')
        with gzip.open(file, "rb") as fp:  # Pickling
            list_of_gene_objects = pickle.load(fp)
        return list_of_gene_objects

    @staticmethod
    def is_ID_in_parent_class(ID):
        '''checks wether ID is a Gene or a Protein_Sequence object attribute'''
        parent_class = True
        if ID in ['ENSG_version', 'ENST', 'ENST_version', 'ENSP', 'ENSP_version', 'refseq_NM','refseq_NM_version', 'refseq_NP','refseq_NP_version',
                  'uniprot_accession', 'uniprot_uniparc', 'uniprot_isoform','UCSC_stable_ID','Uniprot_ID','transcript_name']:
            parent_class = False
        return parent_class


    @staticmethod
    def search_through_database_with_known_ID_Type(list_of_gene_objects, dict_of_IDs):
        '''
        Function that searches trough database with gettatribute()
        :param database_list, list_of_IDs
        :return: dictionary of indexes of each element
        '''
        dict_element_indexes = {}
        for element, ID in dict_of_IDs.items():
            if ID == "aminoacid_sequence":
                dict_element_indexes[element] = 'aminoacid_sequence'
                continue
            found = False
            parent_class = Input_flow.is_ID_in_parent_class(ID)
            for index, gene in enumerate(list_of_gene_objects):
                if found:
                    break
                if parent_class:
                    if ID == "gene_name":
                        if getattr(gene, "ensembl_gene_symbol") == element.upper():
                            dict_element_indexes[element] = index
                            break
                        if getattr(gene, "HGNC_gene_symbol") == element.upper():
                            dict_element_indexes[element] = index
                            break
                        if type(gene.previous_symbols) == list:  # line can be deleted since all these attributes should be lists
                            if element.upper() in getattr(gene, "previous_symbols"):
                                dict_element_indexes[element] = index
                                break
                        if type(gene.alias_symbols) == list:
                            if element.upper() in getattr(gene, "alias_symbols"):
                                dict_element_indexes[element] = index
                                break
                    elif ID =="refseq_gene_ID":
                        if getattr(gene, ID) == float(element):
                            dict_element_indexes[element] = index
                            break
                    else:
                        if getattr(gene, ID) == element:
                            dict_element_indexes[element] = index
                            break
                else:
                    if type(gene.protein_sequence_isoform_collection) == list:
                        for protein_sequence in gene.protein_sequence_isoform_collection:
                            if getattr(protein_sequence, ID) == element:
                                dict_element_indexes[element] = index
                                found = True
                                break
            else:
                dict_element_indexes[element] = 'not found'
        return dict_element_indexes


    @staticmethod
    def show_which_elements_were_not_found(input1_IDs):
        '''
        create streamlit notification of which genes were found and which not
        :param input1_IDs:
        :return: message which contains the elements which were not identified
        '''
        number_of_elements = len(input1_IDs)
        nomatch = 0
        list_of_unmatched_elements = []
        for element,index in input1_IDs.items():
            if index == 'not found' or index == "aminoacid_sequence":
                list_of_unmatched_elements.append(element)
                nomatch +=1
        matched_elements = number_of_elements-nomatch
        if nomatch ==0:
                if number_of_elements >1:
                    pass
                    #st.success('All '+ str(number_of_elements)+' elements were successfully identified.')
                else:
                    pass
                    #st.success('Element succesfully identified')
        elif matched_elements==0:
            if number_of_elements==1 and list(input1_IDs.values())[0]=='aminoacid_sequence':
                st.warning("Looks like an amino acid sequence. Click on 'Clear' and 'add 2nd sequence'")
            else:
                st.warning ('No references found in the library')
        else:
            st.warning(str(matched_elements)+'/'+str(number_of_elements)+' genes were successfully found. Unidentified input: '+', '.join(list_of_unmatched_elements))


    @staticmethod
    def show_which_elements_are_not_canonical_and_one_isoform(list_of_gene_objects, nested_dict, dict_of_IDs):
        ''' warn the user that the reference transcript is not specified and automatically chosen
        :param input1_IDs:
        :return:
        '''
        number_of_elements = len(nested_dict)
        notspecified = 0
        list_of_parent_elements = []
        for element, type in dict_of_IDs.items():
            if Input_flow.is_ID_in_parent_class(type) and element in nested_dict.keys():
                list_of_parent_elements.append(element)
                notspecified += 1
        if notspecified!=0:
            if number_of_elements !=1:
                st.info('Note: '+str(notspecified) + '/' + str(number_of_elements) + " reference isoform ID's were not specified in the Input field and are automatically chosen: ("+', '.join(list_of_parent_elements)+')')
            else:
                pass

        one_isoform = 0
        list_of_one_isoform = []
        for element in nested_dict.items():
            gene_index = list(element[1].keys())[0]
            if len(list_of_gene_objects[gene_index].protein_sequence_isoform_collection)==1:
                list_of_one_isoform.append(element[0])
                one_isoform += 1
        if one_isoform!=0:
            if number_of_elements !=1:
                st.info('For '+str(one_isoform) + '/' + str(
                    number_of_elements) + " of the genes, only one protein sequence was found in the library: (" + ', '.join(
                    list_of_one_isoform) + ')')

            else:
                pass


    @staticmethod
    def remove_dict_elements_with_no_gene_object_match(input1_IDs):  # doesnt work, maybe create a whole new dictionary..?, later to be implemented in generate neseted_ dictionary function
        '''
        :param input1_IDs:
        :return: dictionary which the 'not found' elements were removed
        '''
        cleaned_Input1_IDs = dict()
        for element,index in input1_IDs.items():
            if index == "not found" or index =="aminoacid_sequence":
                continue
            else:
                cleaned_Input1_IDs[element] = index
        return cleaned_Input1_IDs

    @staticmethod
    def pick_index_of_canonical_sequence(list_of_gene_objects, index):
        '''
        returns the index of the protein_isofrom of a gene with the longest AA sequence
        :param list_of_gene_objects:
        :return: index of longest protein_isoform object
        '''
        list_of_AA_sequences = [len(protein_isoform.protein_sequence) for protein_isoform in
                                list_of_gene_objects[index].protein_sequence_isoform_collection if
                                protein_isoform.protein_sequence != None]
        index_of_longest_AA = list_of_AA_sequences.index(max(list_of_AA_sequences))
        return index_of_longest_AA


    @staticmethod
    def generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs, dict_element_indexes,
                                                                        list_of_gene_objects):
        '''
        function that returns a dictionary in which each ID has the canonical ID as a value. If user specified ID, then key=value, if not, the canonical sequence must be extracted from the gene object
        :param ID_dictionary:
        :param list_of_gene_objects:
        :return: dictionary of canonical ID's used as a default reference
        '''



        def find_index_of_reference_transcript(element):
            for index_sequence, protein_object in enumerate(list_of_gene_objects[index].protein_sequence_isoform_collection):
                if getattr(protein_object, dict_of_IDs[element]) == element:
                    index_of_reference_sequence = index_sequence
                    break
            return index_of_reference_sequence

        for element, index in dict_element_indexes.items():
            if Input_flow.is_ID_in_parent_class(dict_of_IDs[element]):
                dict_element_indexes[element] = dict({index:Input_flow.pick_index_of_canonical_sequence(list_of_gene_objects,index)})
            else:
                dict_element_indexes[element] = dict({index: find_index_of_reference_transcript(element)})

        return dict_element_indexes


    @staticmethod
    def chose_columns(list_of_gene_objects,nested_dict,dict_of_IDs,run_id,parameter_change):
        st.markdown("### Mapped Amino Acid Positions Table")
        generate_table = False
        if  parameter_change and run_id>1:
            st.write('\n')
            generate_table = st.button('Click here to generate Mapping Table with updated parameters!')

        if  run_id==1 or generate_table or not parameter_change:
            Input_flow.show_which_elements_are_not_canonical_and_one_isoform(list_of_gene_objects,nested_dict, dict_of_IDs)
            container = st.beta_container()
            all = st.checkbox("Select all columns")

            if all:
                selected_options = container.multiselect('Select further columns',
                                                         ['Gene name', 'Ensembl Gene ID (ENSG)',
                                                          'Ensembl Transcript ID (ENST)', 'Ensembl Protein ID (ENSP)',
                                                          'Transcript Name',
                                                          'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)',
                                                          'Refseq Protein ID (NP)', 'UCSC Stable ID (uc)',
                                                          'Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID',
                                                          'Uniparc ID',
                                                          'Ensembl Gene ID version (ENSG.Number)',
                                                          'Ensembl Transcript ID version (ENST.Number)',
                                                          'Ensembl Protein ID version (ENSP.Number)',
                                                          'Refseq Transcript ID version (NM.Number)',
                                                          'Refseq Transcript ID version (NP.Number)',
                                                          'HGNC ID (HGNC:Number)'],
                                                         ['Gene name', 'Ensembl Gene ID (ENSG)',
                                                          'Ensembl Transcript ID (ENST)', 'Ensembl Protein ID (ENSP)',
                                                          'Transcript Name',
                                                          'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)',
                                                          'Refseq Protein ID (NP)', 'UCSC Stable ID (uc)',
                                                          'Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID',
                                                          'Uniparc ID',
                                                          'Ensembl Gene ID version (ENSG.Number)',
                                                          'Ensembl Transcript ID version (ENST.Number)',
                                                          'Ensembl Protein ID version (ENSP.Number)',
                                                          'Refseq Transcript ID version (NM.Number)',
                                                          'Refseq Transcript ID version (NP.Number)',
                                                          'HGNC ID (HGNC:Number)'])
            else:
                selected_options = container.multiselect('Select further columns',
                                                         ['Gene name', 'Ensembl Gene ID (ENSG)',
                                                          'Ensembl Transcript ID (ENST)', 'Ensembl Protein ID (ENSP)',
                                                          'Transcript Name',
                                                          'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)',
                                                          'Refseq Protein ID (NP)', 'UCSC Stable ID (uc)',
                                                          'Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID',
                                                          'Uniparc ID',
                                                          'Ensembl Gene ID version (ENSG.Number)',
                                                          'Ensembl Transcript ID version (ENST.Number)',
                                                          'Ensembl Protein ID version (ENSP.Number)',
                                                          'Refseq Transcript ID version (NM.Number)',
                                                          'Refseq Transcript ID version (NP.Number)',
                                                          'HGNC ID (HGNC:Number)'],
                                                         ['Gene name', 'Transcript Name'])

            return selected_options

    @staticmethod
    def generate_download_section(df):
        download, format = st.beta_columns([0.22, 1])
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
            st.markdown(
                Streamlit_community.get_table_download_link(df, 'Click here to download', sep),
                unsafe_allow_html=True)

    @staticmethod
    def generate_random_example(list_of_gene_objects):
        '''for show example button on website'''
        one_example=False
        parent=False
        fulfilled = False
        random_number_of_examples = random.randint(1,2)
        while not fulfilled:
            if random_number_of_examples >1:
                one_example=True
            if one_example:
                parent_child = random.randint(1, 3)
                if parent_child>1:
                    parent=True
                    list_of_attributes = [a for a in dir(Gene.Gene()) if not a.startswith('__') and not a.startswith('list_') and not a.startswith('minimal') and not a.startswith('protein') and not a.startswith('alias') and not a.startswith('previous') and not a.startswith('refseq')]
                    ID_type = random.randint(0,4)
                    gene_index = random.randint(0,len(list_of_gene_objects))
                    example = getattr(list_of_gene_objects[gene_index],list_of_attributes[ID_type])
                    if example!=None:
                        fulfilled = True
                        return example
                else:
                    list_of_attributes = [a for a in dir(Protein_isoform.Protein_isoform('jkfbkjsdbfbkjbbd')) if
                                          not a.startswith('__') and not a.startswith('list_') and not a.startswith(
                                              'protein') and not a.startswith(
                                              'collection')]
                    ID_type = random.randint(0, 15)
                    gene_index = random.randint(0, len(list_of_gene_objects))
                    isoform_number=len(list_of_gene_objects[gene_index].protein_sequence_isoform_collection)
                    example= getattr(list_of_gene_objects[gene_index].protein_sequence_isoform_collection[random.randint(0,isoform_number-1)], list_of_attributes[ID_type])
                    if example != None:
                        fulfilled = True
                        return example
            else:
                return 'EGFR, Q9Y6I3-1, ENSG00000074410, ENSP00000075430.7, HGNC:10728, UPI00022F85F1'


    @staticmethod
    def pop_up_if_one_isoform(list_of_gene_objects,index_gene_object):
        st.info('ℹ️ There is only one protein sequence stored of this gene in the human isoform library.')
        st.markdown(' #### Isoform Sequence:')
        st.write('\n')
        st.text(list_of_gene_objects[index_gene_object].protein_sequence_isoform_collection[0].protein_sequence)
        st.markdown(' #### Associated information in the library:')
        st.write('\n')
        gene_dict = dict(list_of_gene_objects[index_gene_object].__dict__)
        gene_dict.pop('protein_sequence_isoform_collection')
        st.write('Gene:', gene_dict)
        isoform_dict = dict(list_of_gene_objects[index_gene_object].protein_sequence_isoform_collection[0].__dict__)
        isoform_dict.pop('collection_of_exons')
        isoform_dict.pop('protein_sequence')
        st.write('Protein Isoform:', isoform_dict)
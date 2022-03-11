import pickle
from Streamlit_app.Streamlit_community import *
import gzip
import random
from .Gene import *
from IsoAligner_core.Protein_isoform import *

class Input_flow:
    pass

    @staticmethod
    @st.cache(allow_output_mutation=True)
    def import_data_from_github(file):
        with gzip.open(file, "rb") as fp:
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
    def create_dict_for_pairwise_mode(nested_dict, list_of_gene_objects):
        gene_indexes = set([list(isoform_index.keys())[0] for element, isoform_index in nested_dict.items()])
        gene_names = dict()
        for gene_index in gene_indexes:
            gene_names[list_of_gene_objects[gene_index].ensembl_gene_symbol]=None
        for element, isoform_index in nested_dict.items():
            for gene_index in gene_indexes:
                if list(isoform_index.keys())[0]==gene_index:
                    if gene_names[list_of_gene_objects[gene_index].ensembl_gene_symbol]==None:
                       gene_names[list_of_gene_objects[gene_index].ensembl_gene_symbol] = [element]
                    else:
                        gene_names[list_of_gene_objects[gene_index].ensembl_gene_symbol].append(element)
        return gene_names


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
                        if type(gene.previous_symbols) == list:
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
                else:
                    pass
        elif matched_elements==0:
            if number_of_elements==1 and list(input1_IDs.values())[0]=='aminoacid_sequence':
                st.warning("Looks like an amino acid sequence. Click on 'enter 2nd sequence' to add another sequence.")
            else:
                st.warning ('No references found in the library')
        else:
            st.warning(str(matched_elements)+'/'+str(number_of_elements)+' genes were successfully found. Unidentified input: '+', '.join(list_of_unmatched_elements))


    @staticmethod
    def show_identical_elements(nested_dict, list_of_gene_objects):
        gene_dict = Input_flow.create_dict_for_pairwise_mode(nested_dict, list_of_gene_objects)
        for gene, elements in gene_dict.items():
            if nested_dict[elements[0]] == nested_dict[elements[1]]:
                st.info("ℹ️ Isoform ID \""+elements[0]+"\" and \""+elements[1]+"\" are associated with the exact same protein sequence. All corresponding IDs can be found in the Details below.")


    @staticmethod
    def report_mode_of_action(nested_dict):
        '''
        :param nested_dict:
        :return: mode of action as variable: pairwise or one_ID_per_gene and streamlit warnings
        '''
        gene_indexes=[list(isoform_index.keys())[0] for element, isoform_index in nested_dict.items()]
        if len(gene_indexes)==0:
            return 'stop'
        if len(gene_indexes)== len(set(gene_indexes))*2:
            return 'pairwise'
        elif len(gene_indexes)!= len(set(gene_indexes)):
            st.warning('Please enter your input in the format of either one **or** two Isoform IDs per gene.')
            return 'stop'
        else:
            return "one_ID_per_gene"


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
                st.info('Note: '+str(notspecified) + '/' + str(number_of_elements) + " reference isoform IDs were not specified in the Input field and are automatically chosen: ("+', '.join(list_of_parent_elements)+')')
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
    def remove_dict_elements_with_no_gene_object_match(input1_IDs):
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
    def import_attribute_column_dict():
        attribute_column_dict = {'gene_name': 'Gene name', 'ENSG': 'Ensembl Gene ID (ENSG)',
                                 'ENST': 'Ensembl Transcript ID (ENST)',
                                 'ENSP': 'Ensembl Protein ID (ENSP)', 'transcript_name': 'Transcript Name',
                                 'refseq_gene_ID': 'Refseq Gene ID (Number)', 'refseq_NM': 'Refseq Transcript ID (NM)',
                                 'refseq_NP': 'Refseq Protein ID (NP)',
                                 'UCSC_stable_ID': 'UCSC Stable ID (uc)',
                                 'uniprot_name_ID': 'Uniprot Name ID', 'uniprot_accession': 'Uniprot Accession ID',
                                 'uniprot_isoform': 'Uniprot Isoform ID', 'uniprot_uniparc': 'Uniparc ID',
                                 'ENSG_version': 'Ensembl Gene ID version (ENSG.Number)',
                                 'ENST_version': 'Ensembl Transcript ID version (ENST.Number)',
                                 'ENSP_version': 'Ensembl Protein ID version (ENSP.Number)',
                                 'refseq_NM_version': 'Refseq Transcript ID version (NM.Number)',
                                 'refseq_NP_version': 'Refseq Transcript ID version (NP.Number)',
                                 'HGNC': 'HGNC ID (HGNC:Number)'}
        return attribute_column_dict


    @staticmethod
    def inlcude_id_type_of_user_input_in_df(dict_of_IDs):
        attribute_column_dict = Input_flow.import_attribute_column_dict()
        ids_in_user_input = [ attribute_column_dict[id_type] for id_type in list(dict_of_IDs.values()) if not id_type in ['gene_name', 'transcript_name','aminoacid_sequence']]
        final_id_list = ['Gene name'] + ids_in_user_input + ['Transcript Name']
        return final_id_list

    @staticmethod
    def chose_columns(list_of_gene_objects,nested_dict,dict_of_IDs,run_id,parameter_change):
        st.markdown("### Mapped Amino Acid Positions Table")
        generate_table = False
        if  parameter_change and run_id>1:
            st.write('\n')
            generate_table = st.button('Generate Mapping Table with updated Needleman-Wunsch parameters')
            st.info('⚠️ Manual gene-specific minimal exon length value **changes** are not taken into account for the computing of the mapping table when investigating isoforms from multiple genes simultaneously. See Manual & About for more details.')

        if run_id==1 or generate_table or not parameter_change:
            Input_flow.show_which_elements_are_not_canonical_and_one_isoform(list_of_gene_objects,nested_dict, dict_of_IDs)
            container = st.container()
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
                                                         Input_flow.inlcude_id_type_of_user_input_in_df(dict_of_IDs))
            return selected_options

    @staticmethod
    def generate_download_section(df, no_column=False):
        st.write('\n')
        if not no_column:
            download, format = st.columns([0.25, 1])
            with download:
                sep_choice = st.radio('Choose file format:', ['tsv', 'csv'])
                if sep_choice == "tsv":
                    sep = '\t'
                else:
                    sep = ','
            with format:
                st.write('\n')
                st.write('\n')
                Streamlit_community.download_button(df, 'AA_mapping_table.'+sep_choice, '⇩ Download Table', df, sep)
        else:
            sep_choice = st.radio('Choose file format:', ['tsv', 'csv'])
            if sep_choice == "tsv":
                sep = '\t'
            else:
                sep = ','
            Streamlit_community.download_button(df, 'AA_mapping_table.' + sep_choice, '⇩ Download Table', df, sep)


    @staticmethod
    def generate_random_example(list_of_gene_objects):
        '''for show example button on website'''
        one_example=False
        parent=False
        fulfilled = False
        random_number_of_examples = random.randint(1,5)
        while not fulfilled:
            if random_number_of_examples >1:
                one_example=True
            if one_example:
                parent_child = random.randint(1, 5)
                if parent_child<2:
                    parent=True
                    list_of_attributes = [a for a in dir(Gene()) if not a.startswith('__') and not a.startswith('list_') and not a.startswith('minimal') and not a.startswith('protein') and not a.startswith('alias') and not a.startswith('previous') and not a.startswith('refseq')]
                    ID_type = random.randint(0,4)
                    gene_index = random.randint(0,len(list_of_gene_objects))
                    if len(list_of_gene_objects[gene_index].protein_sequence_isoform_collection)==1:
                        return 'RPS6KA4-201'
                    example = getattr(list_of_gene_objects[gene_index],list_of_attributes[ID_type])
                    if example!=None:
                        fulfilled = True
                        return example
                    else:
                        return 'RPS6KA4-201'
                else:
                    list_of_attributes = [a for a in dir(Protein_isoform('jkfbkjsdbfbkjbbd')) if
                                          not a.startswith('__') and not a.startswith('list_') and not a.startswith(
                                              'protein') and not a.startswith(
                                              'collection')]
                    ID_type = random.randint(0, 15)
                    gene_index = random.randint(0, len(list_of_gene_objects))
                    isoform_number=len(list_of_gene_objects[gene_index].protein_sequence_isoform_collection)
                    if isoform_number ==1:
                        return 'RPS6KA4-201'
                    example= getattr(list_of_gene_objects[gene_index].protein_sequence_isoform_collection[random.randint(0,isoform_number-1)], list_of_attributes[ID_type])
                    if example != None:
                        fulfilled = True
                        return example
                    else:
                        return 'RPS6KA4-201'
            else:
                return 'RPS6KA4, KRAS, Q9Y6I3-1, ENSG00000074410, ENSP00000075430.7, HGNC:10728, UPI00022F85F1'


    @staticmethod
    def pop_up_isoform_info(list_of_gene_objects, index_gene_object,one_isoform, index_isoform=0):
        if one_isoform:
            st.info('ℹ️ There is only one protein sequence stored of this gene in the human isoform library.')
        st.markdown(' #### Isoform Sequence:')
        st.write('\n')
        AA_seq =  list_of_gene_objects[index_gene_object].protein_sequence_isoform_collection[0].protein_sequence
        AA_desc=Input_flow.generate_position_description(
           AA_seq,
            whitespace="", orentiation='bottom', bottom_length="")
        length_str_top = "length: " + str([1 if item != "-" else 0 for item in AA_seq].count(1)) + " AA"
        st.text(length_str_top)
        st.text(''.join(AA_seq) + "\n" + ''.join(AA_desc))
        st.markdown(' #### Associated information in the library:')
        st.write('\n')
        gene_dict = dict(list_of_gene_objects[index_gene_object].__dict__)
        gene_dict.pop('protein_sequence_isoform_collection')
        if gene_dict['previous_symbols']!= list:
            gene_dict.pop('previous_symbols')
        if gene_dict['alias_symbols'] != list:
            gene_dict.pop('alias_symbols')
        list_of_to_be_popped_gene = [attribute for attribute, value in gene_dict.items() if value ==None]
        for attribute in list_of_to_be_popped_gene:
            gene_dict.pop(attribute)
        st.write('Gene Attributes:', gene_dict)
        isoform_dict = dict(list_of_gene_objects[index_gene_object].protein_sequence_isoform_collection[index_isoform].__dict__)
        isoform_dict.pop('gene_name')
        isoform_dict.pop('ENSG')
        isoform_dict.pop('ENSG_version')
        isoform_dict.pop('collection_of_exons')
        isoform_dict.pop('protein_sequence')
        list_of_to_be_popped = [attribute for attribute, value in isoform_dict.items() if value ==None]
        for attribute in list_of_to_be_popped:
            isoform_dict.pop(attribute)
        st.write('Isoform Attributes:', isoform_dict)



#copied from visualise alignment module because import did not work...

    @staticmethod
    def generate_position_description(sequence_list, whitespace, orentiation="top", bottom_length=None):
        symbols= ["|" if int(index+1) % 10==0 or index==0 else "." for index,element in enumerate(sequence_list)]
        positions = []
        one_digit = False
        two_digit = False
        three_digit = False
        four_digit = False
        AA_numbers = Input_flow.compute_AA_number_for_positions(sequence_list)
        iterator = 0
        for index, symbol in enumerate(symbols):
            if symbol=="|":
                try: #needs fixing, problem at the end of the sequence
                    positions.append(str(AA_numbers[iterator]))
                except:
                    positions.append(str(''))
                iterator +=1
                try:
                    if AA_numbers[iterator-1] <10:
                        one_digit = True
                    if AA_numbers[iterator-1] >= 10:
                        two_digit = True
                    if AA_numbers[iterator-1] >=100:
                        three_digit=True
                    if AA_numbers[iterator-1] >=1000:
                        four_digit=True
                except:
                    one_digit=True
            else:
                if four_digit:
                    four_digit=False
                    continue
                if three_digit:
                    three_digit=False
                    continue
                if two_digit:
                   two_digit=False
                   continue
                if one_digit:
                    positions.append(" ")
                    one_digit= False
                else:
                    positions.append(" ")



        if orentiation=="top":
            construct = ''.join(positions) + '\n' + whitespace + ''.join(symbols)
        elif orentiation=="bottom":
            construct = whitespace + ''.join(symbols) + '\n' + bottom_length + ''.join(positions)

        return construct

    @staticmethod
    def compute_AA_number_for_positions(sequence_list):
        AA_numbers = []
        counter = 0
        for index, element in enumerate(sequence_list):
            if index==0 and element!="-":
                counter +=1
            if index % 10==0:
                if index==0:
                    if not element =="-":
                        AA_numbers.append(counter)
                    else:
                        AA_numbers.append(0)
                else:
                    AA_numbers.append(counter)
            if element != "-" and index != 0:
                counter += 1
        #return Visualise_Alignment.clean_up_AA_numbers(AA_numbers)
        return (AA_numbers)

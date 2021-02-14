import re
import pickle
import streamlit as st
import urllib

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
    @st.cache(allow_output_mutation=True)
    def import_data_from_github(file):
        '''import reference file (database), a pickle file generated in the database.py file'''
        with open(file, "rb") as fp:  # Pickling
            list_of_gene_objects = pickle.load(fp)
        return list_of_gene_objects

    @staticmethod
    def is_ID_in_parent_class(ID):
        '''checks wether ID is a Gene or a Protein_Sequence object attribute'''
        parent_class = True
        if ID in ['ENSG_version', 'ENST', 'ENST_version', 'ENSP', 'ENSP_version', 'refseq_rna', 'refseq_protein',
                  'uniprot_accession', 'uniprot_uniparc', 'uniprot_isoform']:  # list must be completed
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
            found = False
            parent_class = Input_flow.is_ID_in_parent_class(ID)
            for index, gene in enumerate(list_of_gene_objects):
                if found:
                    break
                if parent_class:
                    if ID == "gene_name":
                        if getattr(gene, "ensembl_gene_symbol") == element:
                            dict_element_indexes[element] = index
                            break
                        if getattr(gene, "HGNC_gene_symbol") == element:
                            dict_element_indexes[element] = index
                            break
                        if type(
                                gene.previous_symbols) == list:  # line can be deleted since all these attributes should be lists
                            if element in getattr(gene, "previous_symbols"):
                                dict_element_indexes[element] = index
                                break
                        if type(gene.alias_symbols) == list:
                            if element in getattr(gene, "alias_symbols"):
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
            if index == 'not found':
                list_of_unmatched_elements.append(element)
                nomatch +=1
        matched_elements = number_of_elements-nomatch
        if nomatch ==0:
                st.success('All elements were successfully identified.')
        elif matched_elements==0:
            st.warning ('No references found in the database')
        else:
            st.info(str(matched_elements)+'/'+str(number_of_elements)+' elements were successfully found. Unidentified elements: '+', '.join(list_of_unmatched_elements))


    @staticmethod
    def remove_dict_elements_with_no_gene_object_match(input1_IDs):  # doesnt work, maybe create a whole new dictionary..?, later to be implemented in generate neseted_ dictionary function
        '''
        :param input1_IDs:
        :return: dictionary which the 'not found' elements were removed
        '''
        cleaned_Input1_IDs = dict()
        for element,index in input1_IDs.items():
            if index == "not found":
                continue
            else:
                cleaned_Input1_IDs[element] = index
        return cleaned_Input1_IDs


    @staticmethod
    def generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs, dict_element_indexes,
                                                                        list_of_gene_objects):
        '''
        function that returns a dictionary in which each ID has the canonical ID as a value. If user specified ID, then key=value, if not, the canonical sequence must be extracted from the gene object
        :param ID_dictionary:
        :param list_of_gene_objects:
        :return: dictionary of canonical ID's used as a default reference
        '''

        def pick_index_of_canonical_sequence(index_of_gene_object):
            '''
            returns the index of the protein_isofrom of a gene with the longest AA sequence
            :param list_of_gene_objects:
            :return: index of longest protein_isoform object
            '''
            list_of_AA_sequences = [protein_isoform.protein_sequence for protein_isoform in list_of_gene_objects[index].protein_sequence_isoform_collection]
            index_of_longest_AA = list_of_AA_sequences.index(max(list_of_AA_sequences))
            return index_of_longest_AA

        def find_index_of_reference_transcript(element):
            for index_sequence, protein_object in enumerate(list_of_gene_objects[index].protein_sequence_isoform_collection):
                if getattr(protein_object, dict_of_IDs[element]) == element:
                    index_of_reference_sequence = index_sequence
                    break
            return index_of_reference_sequence

        for element, index in dict_element_indexes.items():
            if Input_flow.is_ID_in_parent_class(dict_of_IDs[element]):
                dict_element_indexes[element] = dict({index:pick_index_of_canonical_sequence(index)})
            else:
                dict_element_indexes[element] = dict({index: find_index_of_reference_transcript(element)})

        return dict_element_indexes
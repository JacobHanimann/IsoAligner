import re
from Alignment import *

class Input_preparation:
    pass

    @staticmethod
    def split_elements_from_user_input_string(string):
        '''
        Function that separates gene names/ID from each other
        :param string:
        :return: list of elements
        '''
        string = string.upper()
        if "\n" in string:
            string = re.sub(' ', '', string)
            list_of_elements = list(string.split('\n'))
        elif "," in string:
            string = re.sub(' ', '', string)
            list_of_elements = list(string.split(','))
        elif "\t" in string:
            list_of_elements = list(string.split('\t'))
        elif " " in string and not 'HGNC' in string:
            list_of_elements = list(string.split(' '))
        else:
            list_of_elements = [string]

        return list_of_elements


    @staticmethod
    def identify_IDs_from_user_text_input(string):
        '''
        Function that identifies which ID's the user typed in with regex. Returns a dict of ID_types which can be used to search through the database more efficiently
        :param list of elements:
        :return: dict of ID_types
        comment: should be updated for all kinds of ID in database
        '''
        dict_of_IDs = {}
        list_of_elements = Input_preparation.split_elements_from_user_input_string(string)
        for element in list_of_elements:
            # ensembl
            if re.search('ENSG\d+\.\d+', element):
                dict_of_IDs[element] = 'ENSG_version'
                continue
            elif re.search('ENSG\d{11}', element):
                dict_of_IDs[element] = 'ENSG'
                continue
            if re.search('ENST\d+\.\d+', element):
                dict_of_IDs[element] = 'ENST_version'
                continue
            elif re.search('ENST\d{11}', element):
                dict_of_IDs[element] = 'ENST'
                continue
            if re.search('ENSP\d+\.\d+', element):
                dict_of_IDs[element] = 'ENSP_version'
                continue
            elif re.search('ENSP\d{11}', element):
                dict_of_IDs[element] = 'ENSP'
                continue
            # uniprot
            if re.search(
                    '[OPQ][0-9][0-9A-Z]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z,0-9]{2}[0-9]|[A-N,R-Z][0-9][A-Z][A-Z,0-9]{2}[0-9][A-Z][A-Z,0-9]{2}[0-9]',
                    element) and '-' in element:
                dict_of_IDs[element] = 'uniprot_isoform'
                continue
            #ensembl
            if re.search('.*-\d+', element):
                dict_of_IDs[element] = 'transcript_name'
                continue
            # uniprot
            if re.search(
                    '[OPQ][0-9][0-9A-Z]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z,0-9]{2}[0-9]|[A-N,R-Z][0-9][A-Z][A-Z,0-9]{2}[0-9][A-Z][A-Z,0-9]{2}[0-9]',
                    element):
                dict_of_IDs[element] = 'uniprot_accession'
                continue
            if re.search('UPI[0-9A-F]+', element):
                dict_of_IDs[element] = 'uniprot_uniparc'
                continue

            if re.search('_HUMAN', element):
                dict_of_IDs[element] = 'uniprot_name_ID'
                continue
            #HGNC
            if re.search('HGNC:\d+', element):
                dict_of_IDs[element] = 'HGNC'
                continue

            # refseq
            if re.search('NM_\d+\.\d+', element):
                dict_of_IDs[element] = 'refseq_NM_version'
                continue
            elif re.search('NM_\d+', element):
                dict_of_IDs[element] = 'refseq_NM'
                continue
            if re.search('NP_\d+\.\d+', element):
                dict_of_IDs[element] = 'refseq_NP_version'
                continue
            elif re.search('NP_\d+', element):
                dict_of_IDs[element] = 'refseq_NP'
                continue
            #UCSC
            if re.search('uc\d+.*', element):
                dict_of_IDs[element] = 'UCSC_stable_ID'
                continue

            # refseq again
            if re.search('\d+', element) and not re.search('[A-Z]',element):
                dict_of_IDs[element] = 'refseq_gene_ID'
                continue

            if Alignment.extract_only_AA_of_Fasta_file(element)!=None:
                dict_of_IDs[element] = 'aminoacid_sequence'

            else:# if no ID's were found, the string is probably a gene name
                dict_of_IDs[element] = 'gene_name'

        return dict_of_IDs
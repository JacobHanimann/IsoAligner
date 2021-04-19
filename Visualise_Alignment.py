import re
from collections import Counter
import streamlit as st
from Alignment import *

class Visualise_Alignment:
    pass

    @staticmethod
    def create_list_gene_selection(list_of_gene_objects, input1_IDs):
        '''
        returns list of gene names, associated isoform count and user input
        also, looks like a nightmare
        :param list_of_gene_objects:
        :param input1_IDs:
        :return: list for the gene selection box
        '''
        gene_list = \
         [list_of_gene_objects[list(index.keys())[0]].ensembl_gene_symbol +
         ' (' +
         str(len(list_of_gene_objects[list(index.keys())[0]].protein_sequence_isoform_collection)) +
         ' Entries) | '
         + element
          if element != list_of_gene_objects[list(index.keys())[0]].ensembl_gene_symbol else
          list_of_gene_objects[list(index.keys())[0]].ensembl_gene_symbol +
          ' (' +
          str(len(list_of_gene_objects[list(index.keys())[0]].protein_sequence_isoform_collection)) +
          ' Entries)'
         for element, index in input1_IDs.items()]

        return gene_list


    @staticmethod
    def clean_chosen_gene(chosen_gene):
        '''clean chosen gene from 'select Gene' string for further use back-end'''
        if '|' in chosen_gene:
            chosen_gene_cleaned = re.split(' \| ', chosen_gene)[1]
        elif "(" in chosen_gene:
            chosen_gene_cleaned = re.split(' \(', chosen_gene)[0]
        else: #one ID
            chosen_gene_cleaned = chosen_gene
        return chosen_gene_cleaned


    @staticmethod
    def fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects, dict_element_indexes, chosen_gene):
        '''
        function that fetches isoform IDs in form of a list with the canonical ID as the first element
        :param  list_of_gene_objects, index_of_gene, optional:ID_type:
        :return: list
        '''
        chosen_gene_cleaned = Visualise_Alignment.clean_chosen_gene(chosen_gene)

        index_of_gene_object = list(dict_element_indexes[chosen_gene_cleaned].keys())[0]
        index_of_canonical_transcript = dict_element_indexes[chosen_gene_cleaned][index_of_gene_object]

        list_of_transcripts = []
        index_count = 0
        for sequence in list_of_gene_objects[index_of_gene_object].protein_sequence_isoform_collection:
            canonical=False
            if index_count==index_of_canonical_transcript:
                canonical = True
            if getattr(sequence, 'transcript_name') !=None:
                if not canonical:
                    list_of_transcripts.append((sequence.transcript_name,index_count))
                else: canonical_element =[(sequence.transcript_name,index_count)]
            elif getattr(sequence, 'ENSP') !=None:
                 if not canonical:
                    list_of_transcripts.append((sequence.ENSP,index_count))
                 else:
                     canonical_element = [(sequence.ENSP,index_count)]
            elif getattr(sequence, 'uniprot_isoform') != None:
                if not canonical:
                    list_of_transcripts.append((sequence.uniprot_isoform,index_count))
                else:
                    canonical_element = [(sequence.uniprot_isoform,index_count)]
            elif getattr(sequence, 'refseq_NP') != None:
                if not canonical:
                    list_of_transcripts.append((sequence.refseq_NP,index_count)) #additional IDs have to be added
                else:
                    canonical_element = [(sequence.refseq_NP, index_count)]
            else:
                list_of_attributes = [a for a in dir(sequence) if not a.startswith('__') and not a.startswith('gene_name') and not a.startswith("protein_sequence")]
                for attribute in list_of_attributes:
                    if getattr(sequence,attribute)!=None:
                        if not canonical:
                            list_of_transcripts.append((getattr(sequence,attribute), index_count))  # additional IDs have to be added
                        else:
                            canonical_element = [(getattr(sequence,attribute), index_count)]
                        break
            index_count +=1

        final_transcript_list = canonical_element + list_of_transcripts
        return final_transcript_list, index_of_gene_object


    @staticmethod
    def calculate_percentage_of_mapped_positions(isoform_check, reference_protein_sequence, isoform_protein_sequence):
        '''
        function that returns the percentage of correctly mapped AA per sequence
        :return: percentage of mapped positions
        '''
        positions_total_reference = len(reference_protein_sequence)
        positions_total_isoform = len(isoform_protein_sequence)
        counter_dict = Counter(isoform_check)
        matches = counter_dict['correct']
        percentage_reference = matches / positions_total_reference
        percentage_isoform = matches / positions_total_isoform
        return percentage_reference, percentage_isoform


    @staticmethod
    @st.cache(allow_output_mutation=True)
    def visualise_alignment_dynamically(reference_sequence_list, isoform_sequence_list, AA_match_evalutation_list,
                                        percentage_reference, percentage_isoform, sequence1='sequence1',
                                        sequence2='sequence2', ):
        '''Function that returns an alignment of two sequences according to the AA_match_evalutation list generated from
         the check_for_wrong_exon_alignments() function in a visually pleasing fashion
         Input: 3 lists, 2 floats and 2 strings
         Output: 1 String (formatted with whitespace and newline character)'''
        correct_match_character = "|"
        wrong_match_character = "x"
        name_percentage_reference_string = sequence1 +":(" + str(round(100 * percentage_reference, 1)) + "%) "
        name_percentage_isoform_string = sequence2 +":(" + str(round(100 * percentage_isoform, 1)) + "%) "
        #fixing whitespace length
        if len(name_percentage_reference_string) >= len(name_percentage_isoform_string):
            whitespace = len(name_percentage_reference_string) * " "
            name_percentage_isoform_string = name_percentage_isoform_string + (len(name_percentage_reference_string)-len(name_percentage_isoform_string)) * " "
        else:
            whitespace = len(name_percentage_isoform_string) * " "
            name_percentage_reference_string = name_percentage_reference_string + (len(name_percentage_isoform_string)-len(name_percentage_reference_string)) * " "
        alignment_character_list = [" " if score == "gap" else correct_match_character if score == "correct" else wrong_match_character for score in AA_match_evalutation_list]
        #final string
        output_alignment_string = name_percentage_reference_string + ''.join(
            reference_sequence_list) + '\n' + whitespace + ''.join(
            alignment_character_list) + "\n" + name_percentage_isoform_string + ''.join(isoform_sequence_list)
        return output_alignment_string


    @staticmethod
    def get_index_of_chosen_transcript(chosen_transcript, transcript_index_list):
        '''
        :return: index of transcript (int) in the Gene.collection_of_protein_sequences attribute
        '''
        for transcript in transcript_index_list:
            if chosen_transcript == transcript[0]:
                return transcript[1] #second element of tuple is the index of the transcript in the collection of protein sequences attribute


    @staticmethod
    def fetch_transcript_name_from_selection_of_attributes(list_of_gene_objects,index_of_gene_object,index_of_transcript):
        '''
        checks which transcript names exists and returns hierarchical
        :return: transcript name (string)
        Note: selection of attributes has to be in the same order as in the fetch_Isoform_IDs_of_sequence_collection() function
        '''
        sequence = list_of_gene_objects[index_of_gene_object].protein_sequence_isoform_collection[index_of_transcript]
        if getattr(sequence, 'transcript_name') != None: #has to be complemented with all attributes
            return getattr(sequence, 'transcript_name')
        elif getattr(sequence, 'ENSP') != None:
            return  getattr(sequence, 'ENSP')
        elif getattr(sequence, 'uniprot_isoform') != None:
            return getattr(sequence, 'uniprot_isoform')
        elif getattr(sequence, 'refseq_NP') != None:
            return getattr(sequence, 'refseq_NP')
        else:
            list_of_attributes = [a for a in dir(sequence) if not a.startswith('__') and not a.startswith('gene_name') and not a.startswith("protein_sequence")]
            for attribute in list_of_attributes:
                if getattr(sequence, attribute) != None:
                    return getattr(sequence, attribute)


    @staticmethod
    def display_alignment_for_one_gene_from_database(index_of_reference_transcript, list_of_gene_objects, index_of_gene, match,
                                                     mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA):
        '''
        executes the visualisation of the alignments for a gene with a given reference_transcript
        :param reference_transcript:
        :param list_of_gene_objects:
        :param index_of_gene:
        :return: streamlit write and text commands

        set reference protein sequence
        iterate through others
        create a function for names of sequence1 and sequence2, has to be the same sequence as generating the transcript_index_list for consistency

        '''
        #get reference AA sequence
        reference_protein_sequence = getattr(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection[index_of_reference_transcript], "protein_sequence")

        transcript_number = 1
        sequence1_name = Visualise_Alignment.fetch_transcript_name_from_selection_of_attributes(list_of_gene_objects,index_of_gene,index_of_reference_transcript)
        for index, transcript in enumerate(list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection):
            if index == index_of_reference_transcript: #skip reference AA
                continue

            #compute alignment
            isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = Alignment.map_AA_Needleman_Wunsch_with_exon_check(
                reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,
                gap_extension_penalty, exon_length_AA)[4:7]

            #calculate statistics of alignment output
            percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check,
                                                                                                                    reference_protein_sequence,
                                                                                                        transcript.protein_sequence)
            st.write('Alignment ' + str(transcript_number))
            #extract transcript name of alternative_ID isoform
            sequence2_name = Visualise_Alignment.fetch_transcript_name_from_selection_of_attributes(list_of_gene_objects,index_of_gene,index)
            st.text(Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta, alignment_isoform_fasta,
                                                    isoform_pattern_check, percentage_reference, percentage_isoform,
                                                    sequence1=sequence1_name, sequence2=sequence2_name))
            st.text('\n')
            transcript_number += 1
import re
from collections import Counter
import streamlit as st
from Alignment import *

class Visualise_Alignment:
    pass

    @staticmethod
    def fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects, dict_element_indexes, chosen_gene, ID="ENSP"):
        '''
        function that fetches isoform IDs in form of a list with the canonical ID as the first element
        :param  list_of_gene_objects, index_of_gene, optional:ID_type:
        :return: list
        '''
        if len(dict_element_indexes) > 1:
            chosen_gene_cleaned = re.split(' \(', chosen_gene)[0]
        else:
            chosen_gene_cleaned = chosen_gene
        index_of_gene_object = list(dict_element_indexes[chosen_gene_cleaned].keys())[0]
        index_of_canonical_transcript = dict_element_indexes[chosen_gene_cleaned][index_of_gene_object]
        list_of_transcripts = [getattr(sequence, ID) for sequence in
                               list_of_gene_objects[index_of_gene_object].protein_sequence_isoform_collection if
                               getattr(sequence, ID) != getattr(
                                   list_of_gene_objects[index_of_gene_object].protein_sequence_isoform_collection[
                                       index_of_canonical_transcript], ID)]
        list_of_transcripts = [getattr(list_of_gene_objects[index_of_gene_object].protein_sequence_isoform_collection[
                                           index_of_canonical_transcript], ID)] + list_of_transcripts
        return list_of_transcripts


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
    def visualise_alignment_dynamically(reference_sequence_list, isoform_sequence_list, AA_match_evalutation_list,
                                        percentage_reference, percentage_isoform, sequence1='sequence1',
                                        sequence2='sequence2', ):
        '''Function that returns an alignment of two sequences according to the AA_match_evalutation list generated from
         the check_for_wrong_exon_alignments() function in a visually pleasing fashion
         Input: 3 lists, 2 floats and 2 strings
         Output: 1 String (formatted with whitespace and newline character)'''
        correct_match_character = "|"
        wrong_match_character = "x"
        alignment_character_list = [
            " " if score == "gap" else correct_match_character if score == "correct" else wrong_match_character for
            score in AA_match_evalutation_list]
        whitespace = len(sequence1) * " " + "         "
        output_alignment_string = sequence1 + ":(" + str(round(100 * percentage_reference, 1)) + '%) ' + ''.join(
            reference_sequence_list) + '\n' + whitespace + ''.join(
            alignment_character_list) + "\n" + sequence2 + ":(" + str(
            round(100 * percentage_isoform, 1)) + '%) ' + ''.join(isoform_sequence_list)
        return output_alignment_string


    @staticmethod
    def display_alignment_for_one_gene_from_database(reference_transcript, list_of_gene_objects, index_of_gene, match,
                                                     mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA,
                                                     ID_type='ENSP'):
        '''
        executes the visualisation of the alignments for a gene with a given reference_transcript
        :param reference_transcript:
        :param list_of_gene_objects:
        :param index_of_gene:
        :return: streamlit write and text commands
        '''
        for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
            if getattr(transcript, ID_type) == reference_transcript:
                reference_protein_sequence = getattr(transcript, "Protein_sequence")
                break
        transcript_number = 1
        for transcript in list_of_gene_objects[index_of_gene].protein_sequence_isoform_collection:
            if getattr(transcript, ID_type) == reference_transcript:
                continue
            isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = Alignment.map_AA_Needleman_Wunsch_with_exon_check(
                reference_protein_sequence, transcript.Protein_sequence, match, mismatch, open_gap_penalty,
                gap_extension_penalty, exon_length_AA)[4:7]
            percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check,
                                                                                                reference_protein_sequence,
                                                                                                transcript.Protein_sequence)
            st.write('Alignment ' + str(transcript_number))
            st.text(Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta, alignment_isoform_fasta,
                                                    isoform_pattern_check, percentage_reference, percentage_isoform,
                                                    sequence1=reference_transcript, sequence2=transcript.ENST))
            st.text('\n')
            transcript_number += 1
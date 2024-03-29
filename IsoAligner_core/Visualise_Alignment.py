from collections import Counter
from IsoAligner_core.Alignment import *
from IsoAligner_core.Input_flow import *

class Visualise_Alignment:
    pass

    @staticmethod
    def create_list_gene_selection(list_of_gene_objects, input1_IDs, pairwise=False):
        '''
        returns list of gene names, associated isoform count and user input
        also, looks like a nightmare
        :param list_of_gene_objects:
        :param input1_IDs:
        :return: list for the gene selection box
        '''
        if not pairwise:
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
        else:
            gene_indexes = [list(isoform_index.keys())[0] for element, isoform_index in input1_IDs.items()]
            return [list_of_gene_objects[gene_index].ensembl_gene_symbol+' (2/'+str(len(list_of_gene_objects[gene_index].protein_sequence_isoform_collection))+" Entries)" for gene_index in set(gene_indexes)]


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
    def fetch_Isoform_IDs_of_sequence_collection(list_of_gene_objects, dict_element_indexes, chosen_gene,dict_of_ids):
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
                    list_of_transcripts.append((sequence.refseq_NP,index_count))
                else:
                    canonical_element = [(sequence.refseq_NP, index_count)]
            else:
                list_of_attributes = [a for a in dir(sequence) if not a.startswith('__') and not a.startswith('gene_name') and not a.startswith("protein_sequence")]
                for attribute in list_of_attributes:
                    if getattr(sequence,attribute)!=None:
                        if not canonical:
                            list_of_transcripts.append((getattr(sequence,attribute), index_count))
                        else:
                            canonical_element = [(getattr(sequence,attribute), index_count)]
                        break
            index_count +=1
        if not Input_flow.is_ID_in_parent_class(dict_of_ids[chosen_gene_cleaned]):
            if list(canonical_element)[0][0] != chosen_gene_cleaned:
               canonical_element = [(list(canonical_element)[0][0]+' | '+chosen_gene_cleaned,list(canonical_element)[0][1])]
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
        mismatch_character = 'X'
        name_percentage_reference_string = sequence1 +":(" + str(round(100 * percentage_reference, 1)) + "%) "
        name_percentage_isoform_string = sequence2 +":(" + str(round(100 * percentage_isoform, 1)) + "%) "
        #fixing whitespace length
        if len(name_percentage_reference_string) >= len(name_percentage_isoform_string):
            whitespace = len(name_percentage_reference_string) * " "
            name_percentage_isoform_string = name_percentage_isoform_string + (len(name_percentage_reference_string)-len(name_percentage_isoform_string)) * " "
        else:
            whitespace = len(name_percentage_isoform_string) * " "
            name_percentage_reference_string = name_percentage_reference_string + (len(name_percentage_isoform_string)-len(name_percentage_reference_string)) * " "

        alignment_character_list = [" " if score == "gap" else correct_match_character if score == "correct" else wrong_match_character if score=='wrong' else mismatch_character for score in AA_match_evalutation_list]
        #final string
        length_str_top = "length:" + str([ 1 if item !="-" else 0 for item in reference_sequence_list].count(1))
        length_str_bottom = "length:" + str([ 1 if item !="-" else 0 for item in isoform_sequence_list].count(1))
        adj_whitespace_top = (len(whitespace)-len(length_str_top)) * " "
        adj_whitespace_bottom = (len(whitespace) - len(length_str_bottom)) * " "
        length_str_bottom_fix = length_str_bottom + adj_whitespace_bottom
        reference_sequence_positions = Visualise_Alignment.generate_position_description(reference_sequence_list,
                                                                                         whitespace, orentiation="top")
        isoform_sequence_positions = Visualise_Alignment.generate_position_description(isoform_sequence_list,
                                                                                       whitespace, orentiation="bottom", bottom_length=length_str_bottom_fix)
        output_alignment_string = length_str_top + adj_whitespace_top + ''.join(
            reference_sequence_positions) + '\n' + name_percentage_reference_string  +''.join(
            reference_sequence_list) + '\n' + whitespace + ''.join(
            alignment_character_list) + "\n" + name_percentage_isoform_string + ''.join(isoform_sequence_list) + '\n' + ''.join(isoform_sequence_positions)
        return output_alignment_string


    @staticmethod
    def generate_position_description(sequence_list, whitespace, orentiation="top", bottom_length=None):
        symbols= ["|" if int(index+1) % 10==0 or index==0 else "." for index,element in enumerate(sequence_list)]
        positions = []
        one_digit = False
        two_digit = False
        three_digit = False
        four_digit = False
        AA_numbers = Visualise_Alignment.compute_AA_number_for_positions(sequence_list)
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

    @staticmethod
    def clean_up_AA_numbers(AA_numbers):
        "function above has to be adjusted to integrate this function"
        cleaned = []
        for index, number in enumerate(AA_numbers):
            if index ==0:
                cleaned.append(number)
            elif number != AA_numbers[index-1]:
                cleaned.append(number)
            else:
                cleaned.append("")
        return cleaned


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
        if getattr(sequence, 'transcript_name') != None:
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
                                                     mismatch, open_gap_penalty, gap_extension_penalty, exon_length_AA, pairwise=None, streamlit=False):
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
            stop= False
            if index == index_of_reference_transcript: #skip reference AA
                continue
            if pairwise:
                for element in pairwise:
                    if index != element[1]:
                        stop = True
                    else:
                        stop=False
                        break

            #compute alignment
            if not stop:
                isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = Alignment.map_AA_Needleman_Wunsch_with_exon_check(
                    reference_protein_sequence, transcript.protein_sequence, match, mismatch, open_gap_penalty,
                    gap_extension_penalty, exon_length_AA, streamlit=True)[4:7]

                #calculate statistics of alignment output
                percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check,
                                                                                                                        reference_protein_sequence,
                                                                                                          transcript.protein_sequence)
                st.write("____________")
                st.write('Alignment ' + str(transcript_number))
                #extract transcript name of alternative_ID isoform
                sequence2_name = Visualise_Alignment.fetch_transcript_name_from_selection_of_attributes(list_of_gene_objects,index_of_gene,index)
                st.text(Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta, alignment_isoform_fasta,
                                                        isoform_pattern_check, percentage_reference, percentage_isoform,
                                                        sequence1=sequence1_name, sequence2=sequence2_name))
                st.text('\n')
                transcript_number += 1
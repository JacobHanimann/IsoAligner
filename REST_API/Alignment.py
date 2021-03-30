import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class Alignment:
    pass

    @staticmethod
    def extract_only_AA_of_Fasta_file(fasta_file):
        '''
        using Regex to match pattern of an AA sequence and creating a list in which each AA is an element
        Input type: String
        Output type: String
        '''
        without_newline = re.sub('\n', '', fasta_file)
        without_newline_and_whitespice = re.sub('\s', '', without_newline)
        sequence_of_AA_acronym = '[A-Z]+'
        minimal_length_of_AA_seq = '[A-Z]{7}'
        raw_AA_seq_list = re.findall(sequence_of_AA_acronym + minimal_length_of_AA_seq, without_newline_and_whitespice)
        if len(raw_AA_seq_list) >= 1:
            return raw_AA_seq_list[0]  # string
        else:
            return None


    @staticmethod
    # This function also needs a sanity check
    def check_for_wrong_exon_alignments(ref, isoform,
                                        exon_length_AA):  # Has to be adjusted so the minimal exon length can be varied with tool (more neighbours have to be counted)
        '''
        This function helps to identify falsely aligned elements (distinct exons) when globally aligning isoforms of a gene.
        The Needleman Wunsch algorithm also randomly alignes fractions of different exons but they do not represent the same aminoacid.
        since the optimization of the algorithm is to maximize matches (which makes sense with homologues but not with isoforms)
        :param ref: aligned sequence in form of a list
        :param isoform: aligned sequence in form of a list
        :return: list which categories each elements alignment in 'correct,wrong,gap'
        To do: Function could be written in fewer lines: few things are copied like category gap classification and appending the category..could be all shortened
        '''
        isoform_check = []
        for index in range(0, len(ref)):
            score = 0  # score which determines if the position is fullfills the minimal exon length parameter
            gap = False
            # start of array
            if index <= exon_length_AA:
                if ref[index] != isoform[index]:
                    category = "gap"
                    gap = True
                else:  # same Aminoacid
                    score += 1
                    for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                        if ref[index + sidestep] == isoform[index + sidestep]:
                            score += 1
                        else:
                            break
                if score >= exon_length_AA and gap != True:
                    category = 'correct'
                elif score <= exon_length_AA and gap != True:
                    category = 'wrong'
                isoform_check.append(category)
                continue
            # end of array
            if len(ref) - index <= exon_length_AA:
                if ref[index] != isoform[index]:
                    category = "gap"
                    gap = True
                else:  # same Aminoacid
                    score += 1
                    for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                        if ref[index - sidestep] == isoform[index - sidestep]:
                            score += 1
                        else:
                            break
                if score >= exon_length_AA and gap != True:
                    category = 'correct'
                elif score <= exon_length_AA and gap != True:
                    category = 'wrong'
                isoform_check.append(category)
                continue
            # middle of array, checks for matches both sided
            if ref[index] != isoform[index]:
                category = "gap"
                gap = True
            else:  # same Aminoacid
                score += 1
                for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                    if ref[index + sidestep] == isoform[index + sidestep]:
                        score += 1
                    else:
                        break
                for sidestep in range(1, exon_length_AA):  # check how many neighbours there are
                    if ref[index - sidestep] == isoform[index - sidestep]:
                        score += 1
                    else:
                        break
            if score >= exon_length_AA and gap != True:
                category = 'correct'
            elif score <= exon_length_AA and gap != True:
                category = 'wrong'
            isoform_check.append(category)
        return isoform_check


    @staticmethod
    def map_AA_Needleman_Wunsch_with_exon_check(fasta1, fasta2, match, mismatch, gap_penalty,
                                                gap_extension_penalty, exon_length_AA):
        'maps FMI AA on COSMIC AA and creates list of AAposition and gaps'
        clean_reference_fasta = Alignment.extract_only_AA_of_Fasta_file(fasta1)
        alignments = pairwise2.align.globalms(clean_reference_fasta, Alignment.extract_only_AA_of_Fasta_file(fasta2), match,
                                              mismatch, gap_penalty, gap_extension_penalty, one_alignment_only=True,
                                              penalize_extend_when_opening=True, penalize_end_gaps=False)
        alignment_reference_fasta = list(alignments[0][0])
        alignment_isoform_fasta = list(alignments[0][1])
        isoform_pattern_check = Alignment.check_for_wrong_exon_alignments(alignment_reference_fasta, alignment_isoform_fasta,
                                                                exon_length_AA)
        reference_position_list = []
        isoform_positions_list = []
        aminoacids = []
        position_reference = 1
        position_isoform = 1
        for i in range(0, len(alignment_reference_fasta)):
            if isoform_pattern_check[i] != 'wrong':
                if alignment_reference_fasta[i] == alignment_isoform_fasta[i]:
                    aminoacids.append(alignment_reference_fasta[i])
                    reference_position_list.append(position_reference)
                    isoform_positions_list.append(position_isoform)
                    position_reference += 1
                    position_isoform += 1
                if alignment_reference_fasta[i] == '-' and alignment_isoform_fasta[i] == '-':  # does it even happen?
                    continue
                if alignment_reference_fasta[i] != '-' and alignment_isoform_fasta[i] == '-':
                    position_reference += 1
                if alignment_reference_fasta[i] == '-' and alignment_isoform_fasta[i] != '-':
                    position_isoform += 1
            else:
                position_reference += 1
                position_isoform += 1
        return (format_alignment(*alignments[0], full_sequences=True), aminoacids, reference_position_list,
                isoform_positions_list, isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta)
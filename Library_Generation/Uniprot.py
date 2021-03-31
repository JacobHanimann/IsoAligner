from Extractions_BioIDs import *
from Protein_isoform import *
from Alignment import *
import re


class Uniprot():
    pass


    @staticmethod
    def add_uniprot_fasta_files(file, list_of_gene_objects):
        '''complement library with fasta sequences from uniprot'''

        # prepare file
        with open(file, "r") as f:
            expenses_txt = f.readlines()
            # Put all the lines into a single string
        whole_txt = "".join(expenses_txt)
        splittext = re.split("\n>", whole_txt)

        # organisation
        no_gene_name = 0
        already_in_accession = 0
        accession_in_but_other_sequence = 0
        already_in_uniprot_isoform = 0
        uniprot_isoform_not_same_seq = 0
        new_isoform_for_gene = 0
        no_gene_match_found = 0
        fasta_count = 0
        protein_sequence_already_without_uniprot_ID = 0

        print(len(splittext))

        # iterate

        for fasta in splittext[1:len(splittext)]:

            # organisation
            gene_name_found = True
            fasta_count += 1
            print(fasta_count)
            uniprot_isoform = False

            # extracting information
            try:
                accession = re.split('\|', fasta)[1]
                if "sp" in fasta:
                    uniprot_name_ID = re.findall(".*_HUMAN", re.split('\|', fasta)[2])[0]
                    sp_human_ID = True
                else:
                    sp_human_ID = False
                if "-" in accession:
                    uniprot_isoform = True
            except:
                print('no accession number')
                print(fasta)
            try:
                gene_name = re.findall("GN=[A-Z,0-9]+", fasta)[0][3:]
            except:
                no_gene_name += 1
                gene_name_found = False
            try:
                protein_sequence = Alignment.extract_only_AA_of_Fasta_file(re.split("\n", fasta, maxsplit=1)[1])
            except:
                print('no AA sequence found')

            # if fasta file contains a gene name
            if gene_name_found:

                # find gene object
                gene_identified = False
                for index, gene in enumerate(list_of_gene_objects):
                    if gene_identified:
                        break
                    list_of_gene_names = [gene.ensembl_gene_symbol] + [gene.HGNC_gene_symbol]
                    if type(gene.alias_symbols) == list:
                        list_of_gene_names = list_of_gene_names + gene.alias_symbols
                    if type(gene.previous_symbols) == list:
                        list_of_gene_names = list_of_gene_names + gene.previous_symbols
                    if gene_name in list_of_gene_names:
                        gene_identified = True
                        break

                # if gene name was found in list of gene objects
                if gene_identified:
                    if sp_human_ID == True:
                        gene.uniprot_name_ID = uniprot_name_ID
                    if type(gene.protein_sequence_isoform_collection) == list:
                        found = False
                        for isoform in gene.protein_sequence_isoform_collection:
                            if found:
                                break
                            if isoform.uniprot_accession == accession:
                                if isoform.protein_sequence == protein_sequence:
                                    already_in_accession += 1
                                    isoform.uniprot_isoform = accession + '-1'
                                    found = True
                                else:
                                    accession_in_but_other_sequence += 1
                                    isoform.uniprot_accession = None  # delete (false) attribute of isoform
                                    gene.protein_sequence_isoform_collection.append(
                                        Protein_isoform(protein_sequence, uniprot_accession=accession,
                                                        uniprot_isoform=accession + "-1", gene_name=gene_name))  # add isoform to collection
                                    found = True

                            elif isoform.uniprot_isoform == accession:
                                if isoform.protein_sequence == protein_sequence:
                                    already_in_uniprot_isoform += 1
                                    found = True
                                else:
                                    uniprot_isoform_not_same_seq += 1
                                    isoform.uniprot_isoform = None  # delete (false) attribute of isoform
                                    gene.protein_sequence_isoform_collection.append(
                                        Protein_isoform(protein_sequence,
                                                        uniprot_accession=Get_Bio_ID.get_bio_IDs_with_regex('uniprot_accession',
                                                                                                 accession),
                                                        uniprot_isoform=accession, gene_name=gene_name))
                                    found = True

                            elif isoform.protein_sequence == protein_sequence:
                                found = True
                                protein_sequence_already_without_uniprot_ID += 1
                                if uniprot_isoform:
                                    isoform.uniprot_accession =Get_Bio_ID.get_bio_IDs_with_regex('uniprot_accession', accession)
                                    isoform.uniprot_isoform = accession
                                else:
                                    isoform.uniprot_accession = accession
                                    isoform.uniprot_isoform = accession + '-1'

                        if not found:
                            new_isoform_for_gene += 1
                            if uniprot_isoform:
                                gene.protein_sequence_isoform_collection.append(Protein_isoform(protein_sequence,
                                                                                                uniprot_accession=Get_Bio_ID.get_bio_IDs_with_regex(
                                                                                                    'uniprot_accession',
                                                                                                    accession),
                                                                                                uniprot_isoform=accession,
                                                                                                gene_name=gene_name))
                            else:
                                gene.protein_sequence_isoform_collection.append(
                                    Protein_isoform(protein_sequence, uniprot_accession=accession,
                                                    uniprot_isoform=accession + '-1', gene_name=gene_name))



                # gene name was not found in list of gene objects
                else:
                    no_gene_match_found += 1
                    pass

        print('fasta files with no gene names:', no_gene_name)
        print('already in accession', already_in_accession)
        print('accession but other sequence', accession_in_but_other_sequence)
        print('already in uniprot isoform', already_in_uniprot_isoform)
        print('uniprot_isoform_not_same_seq', uniprot_isoform_not_same_seq)
        print('new_isoform_for_gene', new_isoform_for_gene)
        print('no_gene_match_found', no_gene_match_found)
        print('protein sequence match but no uniprot_name_ID associated: ', protein_sequence_already_without_uniprot_ID)



from Gene import *
from Protein_isoform import *
from Alignment import *
from Extractions_BioIDs import *


class Ensembl():
    pass

    @staticmethod
    def get_ensembl_fasta_sequences_and_IDs_and_create_gene_objects(file):
        '''extract fasta files one by one and add them to the gene objects'''
        with open(file, "r") as f:
            expenses_txt = f.readlines()
        # Put all the lines into a single string
        whole_txt = "".join(expenses_txt)
        splittext = re.split(">", whole_txt)
        fasta_count = 0
        matches = 0
        list_of_gene_objects = []
        for fasta in splittext[1:]:
            fasta_count += 1
            found = False
            gene_name = get_bio_IDs_with_regex('gene_name', fasta)
            # create Protein_isoform object to add to the gene_object
            aa_sequence = Alignment.extract_only_AA_of_Fasta_file(fasta.split('\n', 1)[1])
            if aa_sequence == None:  # sequence shorter than 7 AA long
                continue
            sequence_object = Protein_isoform(aa_sequence, gene_name,
                                              get_bio_IDs_with_regex('ensembl_ensg', fasta),
                                              get_bio_IDs_with_regex('ensembl_ensg_version', fasta),
                                              get_bio_IDs_with_regex('ensembl_enst', fasta),
                                              get_bio_IDs_with_regex('ensembl_enst_version', fasta),
                                              get_bio_IDs_with_regex('ensembl_ensp', fasta),
                                              get_bio_IDs_with_regex('ensembl_ensp_version', fasta),
                                              uniprot_accession=get_bio_IDs_with_regex('uniprot_accession', fasta),
                                              uniprot_uniparc=get_bio_IDs_with_regex('uniprot_uniparc', fasta))
            for gene in list_of_gene_objects:
                if found:
                    break
                if gene.ENSG == sequence_object.ENSG:
                    if type(gene.protein_sequence_isoform_collection) == list:
                        gene.protein_sequence_isoform_collection.append(sequence_object)
                    else:
                        gene.protein_sequence_isoform_collection = [sequence_object]
                    found = True

            if found == False:
                list_of_gene_objects.append(
                    Gene(sequence_object.ENSG, gene_name, protein_sequence_isoform_collection=[sequence_object]))

            print('Fasta files processed: ' + str(fasta_count) + '/' + str(len(splittext)))
        print('Fasta files matched: ' + str(matches))
        return list_of_gene_objects

    


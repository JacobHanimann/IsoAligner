from Gene import *
from Protein_isoform import *
from Alignment import *



class Validate_library():
    pass


    @staticmethod
    def check_if_there_are_AA_seq_duplicates(list_of_gene_objects):
        '''
        check out if there were IDs and Seq that are the same but escaped the match
        :param list_of_gene_objects:
        :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
        '''

        def duplicates(lst, item):
            return [i for i, x in enumerate(lst) if x == item]

        genes_without_AA_seq = 0
        duplicates_number = 0
        genes_without_duplicates = 0
        genes_with_more_than_one_duplicate = 0
        redundant_sequences = 0
        duplicate_genes_dict = dict()
        for index, gene in enumerate(list_of_gene_objects):
            if type(gene.protein_sequence_isoform_collection) == list:
                List = [sequence.protein_sequence for sequence in gene.protein_sequence_isoform_collection]
                duplicates_dict = dict((x, duplicates(List, x)) for x in set(List) if List.count(x) > 1)
                if len(duplicates_dict) != 0:
                    if list(duplicates_dict.keys())[0] != None:
                        duplicate_genes_dict[index] = duplicates_dict
                        duplicates_number += 1
                        for sequence, objects in duplicates_dict.items():
                            redundant_sequences = redundant_sequences + len(objects)
                    if len(duplicates_dict) > 1:
                        genes_with_more_than_one_duplicate += 1
                else:
                    genes_without_duplicates += 1
            else:
                genes_without_AA_seq += 1

        print('number of genes: ', len(list_of_gene_objects))
        print('genes with no AA seq: ', genes_without_AA_seq)
        print('number of genes with AA seq duplicates: ', duplicates_number)
        print('number of genes without AA seq duplicates: ', genes_without_duplicates)
        print('number of redundant AA sequences:', redundant_sequences)
        print('number of genes with more than one AA seq duplicate: ', genes_with_more_than_one_duplicate)
        return duplicate_genes_dict, duplicates_number, genes_without_duplicates, redundant_sequences, genes_with_more_than_one_duplicate


    @staticmethod
    def check_if_there_are_exact_duplicates(list_of_gene_objects):
        '''
        note: function does not work because exon collection attribute is a list which is not hashable
        check out if there were IDs and Seq that are the same but escaped the match
        :param list_of_gene_objects:
        :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
        '''

        def duplicates(lst, item):
            return [i for i, x in enumerate(lst) if x == item]

        genes_without_AA_seq = 0
        duplicates_number = 0
        genes_without_duplicates = 0
        genes_with_more_than_one_duplicate = 0
        duplicate_genes_dict = dict()
        for index, gene in enumerate(list_of_gene_objects):
            if type(gene.protein_sequence_isoform_collection) == list:
                List = [tuple(list(sequence.__dict__.items())) for sequence in gene.protein_sequence_isoform_collection]
                print(dict((x, duplicates(List, x)) for x in set(List) if List.count(x) > 1))
                duplicates_dict = dict((x, duplicates(List, x)) for x in set(List) if List.count(x) > 1)
                if len(duplicates_dict) != 0:
                    if list(duplicates_dict.keys())[0] != None:
                        duplicate_genes_dict[index] = duplicates_dict
                        duplicates_number += 1
                    if len(duplicates_dict) > 1:
                        genes_with_more_than_one_duplicate += 1
                else:
                    genes_without_duplicates += 1
            else:
                genes_without_AA_seq += 1

        print('number of genes: ', len(list_of_gene_objects))
        print('genes with no AA seq: ', genes_without_AA_seq)
        print('number of genes with exact duplicates: ', duplicates_number)
        print('number of genes without exact duplicates: ', genes_without_duplicates)
        print('number of genes with more than one exact duplicate: ', genes_with_more_than_one_duplicate)
        return duplicate_genes_dict

    @staticmethod
    def check_if_there_gene_names_duplicates_over_all_genes(list_of_gene_objects):
        '''
        check out if there were IDs and Seq that are the same but escaped the match
        :param list_of_gene_objects:
        :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
        '''

        def duplicates(lst, item):
            return [i for i, x in enumerate(lst) if x == item]

        list_of_all_names = []
        for index, gene in enumerate(list_of_gene_objects):
            list_of_all_names.append(gene.ensembl_gene_symbol)
        duplicates_dict_gene = dict(
            (x, duplicates(list_of_all_names, x)) for x in set(list_of_all_names) if list_of_all_names.count(x) > 1)

        print(duplicates_dict_gene)
        print('gene object with same name:', len(duplicates_dict_gene))
        return duplicates_dict_gene


    @staticmethod
    def fuse_gene_objects_with_same_name(list_of_gene_objects, duplicates_dict_gene):
        pass


    @staticmethod
    def fuse_attributes_of_duplicated_AA_seq_within_gene_object(list_of_gene_objects, duplicate_genes_dict):
        '''
        function that fuses protein isoform objects if the attributes can complement each other to one big object, otherwise the duplicates will stay separated.
        :param list_of_gene_objects:
        :param duplicate_genes_dict:
        :return: updated list_of_gene_objects
        '''
        reduced_isoform_count = 0
        couldnotmatch = 0
        duplicates_in_total = 0
        for gene, duplicates_dict in duplicate_genes_dict.items():
            tobedeleted = []
            duplicates_in_total = duplicates_in_total + len(duplicates_dict)
            for duplicate_AA in duplicates_dict.items():
                new_object_attributes = Protein_isoform(duplicate_AA[0])
                isoform_dict = dict()
                list_of_attributes = [a for a in dir(new_object_attributes) if not a.startswith('__')]
                different_attributes = False
                for isoform in duplicate_AA[1]:
                    if different_attributes:
                        break
                    isoform = list_of_gene_objects[gene].protein_sequence_isoform_collection[isoform]
                    for attribute in list_of_attributes:
                        if different_attributes:
                            break
                        if getattr(new_object_attributes, attribute) == None:
                            if getattr(isoform, attribute) != None:
                                setattr(new_object_attributes, attribute, getattr(isoform, attribute))
                        else:
                            if getattr(isoform, attribute) != None:
                                if getattr(isoform, attribute) == getattr(new_object_attributes, attribute):
                                    pass  # attributes are the same, protein object can still be fused
                                else:  # stop process, IDs differ from each other
                                    different_attributes = True
                                    couldnotmatch += 1
                if different_attributes == False:
                    tobedeleted.extend(duplicate_AA[1])
                    list_of_gene_objects[gene].protein_sequence_isoform_collection.append(new_object_attributes)
                    reduced_isoform_count += 1
            if tobedeleted:
                for ele in sorted(tobedeleted, reverse=True):
                    del list_of_gene_objects[gene].protein_sequence_isoform_collection[ele]

        print('duplicates in total:', duplicates_in_total)
        print('duplicates that could not be matched:', couldnotmatch)
        print('duplicates that could be matched:', reduced_isoform_count)
        return list_of_gene_objects


    @staticmethod
    def check_if_gene_name_and_prot_seq_are_switched(list_of_gene_objects):
        '''somewhere in the database generation gene name and protein sequence attribute of a protein isoform object are being falsely switched'''
        false_assigned_gene_name_isoform = 0
        for gene in list_of_gene_objects:
            if type(gene.protein_sequence_isoform_collection) == list:
                for isoform in gene.protein_sequence_isoform_collection:
                    if type(isoform.gene_name) == str:
                        if Alignment.extract_only_AA_of_Fasta_file(isoform.gene_name) != None:
                            false_assigned_gene_name_isoform += 1
        print('number of falsely assigned AA seq to gene_name:', false_assigned_gene_name_isoform)


    @staticmethod
    def delete_genes_and_protein_isoforms_with_no_AA_seq(list_of_gene_objects):
        '''
        function that delets empty gene_objects
        :param list_of_gene_obejcts:
        :return: uptaded list of gene objects
        '''
        tobedeletedgene = []
        one_AA_seq = 0
        for index, gene in enumerate(list_of_gene_objects):
            if type(gene.protein_sequence_isoform_collection) != list:
                tobedeletedgene.append(index)
            else:
                if len(gene.protein_sequence_isoform_collection) == 1:
                    one_AA_seq += 1
        if tobedeletedgene:
            for ele in sorted(tobedeletedgene, reverse=True):
                del list_of_gene_objects[ele]

        deleted = 0
        for index_gene, gene in enumerate(list_of_gene_objects):
            tobedeletedisoform = []
            for index, isoform in enumerate(gene.protein_sequence_isoform_collection):
                if isoform.protein_sequence == None:
                    tobedeletedisoform.append(index)
            if tobedeletedisoform:
                for ele in sorted(tobedeletedisoform, reverse=True):
                    del list_of_gene_objects[index_gene].protein_sequence_isoform_collection[ele]
                    deleted += 1

        print('no AA seq Isoforms deleted:', deleted)
        print('genes with just one isoform:', one_AA_seq)
        return list_of_gene_objects
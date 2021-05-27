import pandas as pd

class Biomart_tables():
    pass


    @staticmethod
    def add_Uniprot_Isoform_refseqrna_transcript_name_ID_to_protein_attributes(file, list_of_gene_objects):
        '''
        add IDs to protein isoform object attributes to extend ID library
        :param file: from Biomart
        :return: updated list_of_gene_objects with more attributes
        '''
        df = pd.read_csv(file, sep='\t')
        print('total length:', len(df))
        for index in range(0, len(df)):
            if index % 100 == 0:
                print(100 * round(index / len(df), 2), '%')
            found = False
            uniparc_ID = df.loc[index, 'UniParc ID']
            if type(uniparc_ID) == float:
                continue
            for gene in list_of_gene_objects:
                if found:
                    break
                if type(gene.protein_sequence_isoform_collection) == list:
                    for sequence in gene.protein_sequence_isoform_collection:
                        if found:
                            break
                        if sequence.uniprot_uniparc == uniparc_ID:
                            found = True
                            if type(df.loc[index, 'RefSeq mRNA ID']) != float:
                                sequence.refseq_NM = df.loc[index, 'RefSeq mRNA ID']
                            if type(df.loc[index, 'Transcript name']) != float:
                                sequence.transcript_name = df.loc[index, 'Transcript name']
                            if type(df.loc[index, 'UniProtKB isoform ID']) != float:
                                sequence.uniprot_isoform = df.loc[index, 'UniProtKB isoform ID']
                            break
                else:
                    continue


    @staticmethod
    def add_UCSC_to_protein_attributes(file, list_of_gene_objects):
        '''
        add IDs to protein isoform object attributes to extend ID library
        :param file: from Biomart
        :return: updated list_of_gene_objects with more attributes
        '''
        df = pd.read_csv(file, sep='\t')
        print('total length:', len(df))
        for index in range(0, len(df)):
            if index % 100 == 0:
                print(100 * round(index / len(df), 2), '%')
            found = False
            uniparc_ID = df.loc[index, 'UniParc ID']
            if type(uniparc_ID) == float:
                continue
            for gene in list_of_gene_objects:
                if found:
                    break
                if type(gene.protein_sequence_isoform_collection) == list:
                    for sequence in gene.protein_sequence_isoform_collection:
                        if found:
                            break
                        if sequence.uniprot_uniparc == uniparc_ID:
                            found = True
                            if type(df.loc[index, 'UCSC Stable ID']) != float:
                                sequence.UCSC_stable_ID = df.loc[index, 'UCSC Stable ID']
                            break  # here the uniprotKB gene names ID could be added
                else:
                    continue



    @staticmethod
    def add_refseq_protein_IDs(file, list_of_gene_objects):
        '''add IDs from Biomart file'''
        df = pd.read_csv(file, sep='\t')
        print('total length', len(df))
        for index in range(0, len(df)):
            if index % 100 == 0:
                print(100 * round(index / len(df), 2), '%')
            found = False
            uniparc_ID = df.loc[index, 'UniParc ID']
            if type(uniparc_ID) == float:
                continue
            for gene in list_of_gene_objects:
                if found:
                    break
                if type(gene.protein_sequence_isoform_collection) == list:
                    for sequence in gene.protein_sequence_isoform_collection:
                        if found:
                            break
                        if sequence.uniprot_uniparc == uniparc_ID:
                            found = True
                            if type(df.loc[index, 'RefSeq peptide ID']) != float:
                                sequence.refseq_protein = df.loc[index, 'RefSeq peptide ID']
                                break
                else:
                    continue
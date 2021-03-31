import pandas as pd
from Gene import *


class HGNC():
    pass


    @staticmethod
    def add_HCGN_information_to_gene_objects(file_of_gene_names, list_of_gene_objects):
        '''complement list of gene objects
        input: list of gene objects
        output: list of gene objects with added attribute values
        '''
        df = pd.read_csv(file_of_gene_names, sep='\t')
        print('total length: ', len(df))
        for index in range(0, len(df)):
            print(index)
            # extract data line by line
            HGNC = df.loc[index, 'HGNC']
            HGNC_gene_symbol = df.loc[index, 'approved_symbol']
            previous_symbols = df.loc[index, 'previous_symbols']
            refseq_gene_ID = df.loc[index, 'NCBI Gene ID']
            alias_symbols = df.loc[index, 'alias_symbols']
            uniprot_ID = df.loc[index, 'UniProt ID']  # check if it the same ID as in the Protein_isoform classes

            # transfrom data in correct format
            if type(previous_symbols) != float:  # None values are type float
                if "," in previous_symbols:
                    previous_symbols = previous_symbols.split(', ')
                else:
                    previous_symbols = [
                        previous_symbols]  # either way create a list because it facilitates later search functions
            if type(alias_symbols) != float:
                if "," in alias_symbols:
                    alias_symbols = alias_symbols.split(', ')
                else:
                    alias_symbols = [alias_symbols]

            found = False
            # search for a match
            for gene in list_of_gene_objects:
                if found:
                    break
                if gene.ENSG == df.loc[index, 'Ensembl gene ID']:
                    found = True
                    gene.HGNC = HGNC
                    gene.HGNC_gene_symbol = HGNC_gene_symbol
                    gene.previous_symbols = previous_symbols
                    gene.refseq_gene_ID = refseq_gene_ID
                    gene.alias_symbols = alias_symbols
                    gene.uniprot_name_ID = uniprot_ID  # check if it the same ID as in the Protein_isoform classes

            if found == False:
                list_of_gene_objects.append(
                    Gene(ENSG=df.loc[index, 'Ensembl gene ID'], HGNC=HGNC, HGNC_gene_symbol=HGNC_gene_symbol,
                         previous_symbols=previous_symbols, alias_symbols=alias_symbols, refseq_gene_ID=refseq_gene_ID))

        return list_of_gene_objects


class Gene:
    def __init__(self, ENSG=None, ensembl_gene_symbol=None, refseq_gene_ID=None, HGNC=None, HGNC_gene_symbol=None, previous_symbols=None, alias_symbols=None, protein_sequence_isoform_collection=None, minimal_exon_length=None, uniprot_name_ID=None):
        self.ENSG = ENSG
        self.ensembl_gene_symbol = ensembl_gene_symbol
        self.refseq_gene_ID = refseq_gene_ID
        self.HGNC = HGNC
        self.HGNC_gene_symbol = HGNC_gene_symbol
        self.previous_symbols = previous_symbols
        self.alias_symbols = alias_symbols
        self.protein_sequence_isoform_collection = protein_sequence_isoform_collection
        self.minimal_exon_length= minimal_exon_length
        self.uniprot_name_ID = uniprot_name_ID


    @staticmethod
    def list_of_attributes():
        list_of_attributes_genes = [a for a in dir(Gene()) if not a.startswith('__') and not a.startswith('list_')]
        return list_of_attributes_genes


class Protein_isoform:
    '''objects stored in protein_sequence_isoform_collection attribute of the Gene Class'''
    def __init__(self,  protein_sequence, gene_name=None, ENSG=None, ENSG_version=None, ENST=None, ENST_version=None, ENSP=None, ENSP_version=None,transcript_name=None,
                 refseq_NM=None, refseq_NM_version=None, refseq_NP_version=None, refseq_NP=None, refseq_NC_version=None,refseq_XM_version=None,
                 refseq_XP=None, refseq_XP_version=None, refseq_YP_version=None, refseq_YP=None,
                 uniprot_accession=None, uniprot_uniparc=None, uniprot_isoform=None, uniprot_ID=None, UCSC_stable_ID=None, collection_of_exons=None):
        self.gene_name= gene_name
        self.protein_sequence = protein_sequence
        #ensembl
        self.ENSG = ENSG
        self.ENSG_version = ENSG_version
        self.ENST = ENST
        self.ENST_version = ENST_version
        self.ENSP = ENSP
        self.ENSP_version = ENSP_version
        self.transcript_name = transcript_name
        #refseq
        self.refseq_NM = refseq_NM
        self.refseq_NM_version = refseq_NM_version
        self.refseq_NP = refseq_NP
        self.refseq_NP_version = refseq_NP_version
        #uniprot
        self.uniprot_accession = uniprot_accession
        self.uniprot_uniparc = uniprot_uniparc
        self.uniprot_isoform = uniprot_isoform
        #UCSC
        self.UCSC_stable_ID = UCSC_stable_ID
        #Exon_information
        self.collection_of_exons = collection_of_exons


    @staticmethod
    def list_of_attributes():
        list_of_attributes_isoform = [a for a in dir(Protein_isoform("AJ")) if not a.startswith('__') and not a.startswith('gene_name') and not a.startswith('protein_seq') and not a.startswith('list')]
        return list_of_attributes_isoform
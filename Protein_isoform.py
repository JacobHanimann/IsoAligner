

class Protein_isoform:
    '''objects stored in protein_sequence_isoform_collection attribute of the Gene Class'''
    def __init__(self,  protein_sequence, gene_name=None, ENSG=None, ENSG_version=None, ENST=None, ENST_version=None, ENSP=None, ENSP_version=None,transcript_name=None,
                 refseq_NM=None, refseq_NM_version=None, refseq_NP_version=None, refseq_NP=None, refseq_NC_version=None,refseq_XM_version=None,
                 refseq_XP=None, refseq_XP_version=None, refseq_YP_version=None, refseq_YP=None,
                 uniprot_accession=None, uniprot_uniparc=None, uniprot_isoform=None):
        self.gene_name= gene_name #maybe unnecessary
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
        self.refseq_NC_version = refseq_NC_version
        self.refseq_YP_version = refseq_YP_version
        self.refseq_YP = refseq_YP
        #uniprot
        self.uniprot_accession = uniprot_accession
        self.uniprot_uniparc = uniprot_uniparc
        self.uniprot_isoform = uniprot_isoform


    @staticmethod
    def list_of_attributes():
        list_of_attributes_isoform = [a for a in dir(Protein_isoform("AJ")) if not a.startswith('__') and not a.startswith('gene_name') and not a.startswith('protein_seq') and not a.startswith('list')]
        return list_of_attributes_isoform
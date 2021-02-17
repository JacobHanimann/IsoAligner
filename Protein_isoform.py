

class Protein_isoform:
    def __init__(self, gene_name, protein_sequence, ENSG=None, ENSG_version=None, ENST=None, ENST_version=None, ENSP=None,
                 ENSP_version=None, refseq_NM=None, refseq_NM_version=None, refseq_NP_version=None, refseq_NP=None, uniprot_accession=None, uniprot_uniparc=None, uniprot_isoform=None, transcript_name=None):
        self.gene_name= gene_name #maybe unnecessary
        self.protein_sequence = protein_sequence
        self.ENSG = ENSG
        self.ENSG_version = ENSG_version
        self.ENST = ENST
        self.ENST_version = ENST_version
        self.ENSP = ENSP
        self.ENSP_version = ENSP_version
        self.refseq_NM = refseq_NM
        self.refseq_NM_version = refseq_NM_version
        self.refseq_NP = refseq_NP
        self.refseq_NP_version = refseq_NP_version
        self.uniprot_accession = uniprot_accession
        self.uniprot_uniparc = uniprot_uniparc
        self.uniprot_isoform = uniprot_isoform
        self.transcript_name = transcript_name
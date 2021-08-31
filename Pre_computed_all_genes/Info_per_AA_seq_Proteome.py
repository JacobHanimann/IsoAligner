from IsoAligner_core.Input_flow import *
import pandas as pd

list_of_gene_objects = Input_flow.import_data_from_github('../Human_Isoform_Library/list_of_gene_objects_25th_july.txt.gz')

all_rows = []
for gene in list_of_gene_objects:
    for sequence in gene.protein_sequence_isoform_collection:

        row = [
        gene.ensembl_gene_symbol,
        gene.ENSG,
        sequence.ENSG_version,
        gene.refseq_gene_ID,
        gene.HGNC,
        gene.uniprot_name_ID,
        gene.minimal_exon_length,
        sequence.ENST,
        sequence.ENST_version,
        sequence.ENSP,
        sequence.ENSP_version,
        sequence.transcript_name,
        sequence.refseq_NM,
        sequence.refseq_NM_version,
        sequence.refseq_NP,
        sequence.refseq_NP_version,
        sequence.uniprot_accession,
        sequence.uniprot_isoform,
        sequence.uniprot_uniparc,
        sequence.UCSC_stable_ID,
        sequence.protein_sequence
        ]

        all_rows.append(row)

df = pd.DataFrame(all_rows, columns=['gene_name','ENSG','ENSG_version','refseq_gene_ID',
                                     'HGNC_ID', 'uniprot_name_ID', 'minimal_exon_length',
                                     'ENST', 'ENST_version', 'ENSP', 'ENSP_version',
                                     'transcript_name_ensembl', 'refseq_NM', 'refseq_NM_version',
                                     'refseq_NP', 'refseq_NP_version',
                                     'uniprot_accession', 'uniprot_isoform',
                                     'uniprot_uniparc','UCSC_stable_ID','protein_sequence'], dtype=float)

print(df)

df.to_csv('/Users/jacob/Desktop/Info_per_AA_seq_human_proteome.tsv', sep='\t')
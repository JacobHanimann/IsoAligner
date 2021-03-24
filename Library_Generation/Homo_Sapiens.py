from Ensembl import *
from HGNC import *
from Biomart_tables import *
from Refseq import *
from Uniprot import *

print('Creating list of gene objects with Ensembl Fasta files...')
list_of_gene_objects =Ensembl.get_ensembl_fasta_sequences_and_IDs_and_create_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/ensembl_fasta_IDs_gene_name.txt')

print('Adding HGNC gene symbols to gene attributes...')
HGNC.add_HCGN_information_to_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/HGNC_protein_coding_ensembl.txt',list_of_gene_objects)

print('Adding IDs from Biomart...')
Biomart_tables.add_Uniprot_Isoform_refseqrna_transcript_name_ID_to_protein_attributes('/Users/jacob/Desktop/Isoform Mapper Webtool/NM_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
Biomart_tables.add_refseq_protein_IDs('/Users/jacob/Desktop/Isoform Mapper Webtool/NP_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
Biomart_tables.add_UCSC_to_protein_attributes('/Users/jacob/Desktop/Isoform Mapper Webtool/UCSC_Uniprot_ID_uniparc.txt',list_of_gene_objects)

print('Adding Fasta files from Refseq...')
Refseq.add_refseq_fasta_sequences('/Users/jacob/Desktop/Isoform Mapper Webtool/refseq_fasta_and_info/GCF_000001405.39_GRCh38.p13_protein.gpff',list_of_gene_objects)

print('Adding Fasta files from Uniprot...')
Uniprot.add_uniprot_fasta_files('/Users/jacob/Desktop/Isoform Mapper Webtool/refseq_fasta_and_info/GCF_000001405.39_GRCh38.p13_protein.gpff',list_of_gene_objects)


print('Validating and Correcting Library...')



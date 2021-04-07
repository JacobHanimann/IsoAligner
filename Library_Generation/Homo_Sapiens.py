from Ensembl import *
from HGNC import *
from Biomart_tables import *
from Refseq import *
from Uniprot import *
from Validation_of_library import *
from median_exon_length import *

date = '1st_april'
#
#print('Creating list of gene objects with Ensembl Fasta files...')
#list_of_gene_objects =Ensembl.get_ensembl_fasta_sequences_and_IDs_and_create_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/ensembl_fasta_IDs_gene_name.txt')
#
#print('Pickling list of gene objects and saving file...')
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_first.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
#print('Adding HGNC gene symbols to gene attributes...')
#HGNC.add_HCGN_information_to_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/HGNC_protein_coding_ensembl.txt',list_of_gene_objects)
#
#print('Pickling list of gene objects and saving file...')
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_second.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
#print('Adding IDs from Biomart...')
#Biomart_tables.add_UCSC_to_protein_attributes('/Users/jacob/Desktop/Isoform Mapper Webtool/UCSC_IDs.txt',list_of_gene_objects)
#Biomart_tables.add_Uniprot_Isoform_refseqrna_transcript_name_ID_to_protein_attributes('/Users/jacob/Desktop/Isoform Mapper Webtool/NM_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
#Biomart_tables.add_refseq_protein_IDs('/Users/jacob/Desktop/Isoform Mapper Webtool/NP_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
#
#print('Pickling list of gene objects and saving file...')
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_third.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
#print('Adding Fasta files from Refseq...')
#Refseq.add_refseq_fasta_sequences('/Users/jacob/Desktop/Isoform Mapper Webtool/refseq_fasta_and_info/GCF_000001405.39_GRCh38.p13_protein.gpff',list_of_gene_objects)
#
#print('Pickling list of gene objects and saving file...')
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_fourth.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
#print('Adding Fasta files from Uniprot...')
#Uniprot.add_uniprot_fasta_files('/Users/jacob/Desktop/Isoform Mapper Webtool/uniprot_downloads/uniprot-proteome_UP000005640.fasta',list_of_gene_objects)
#
#print('Pickling list of gene objects and saving file...')
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_fifth.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)

with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_fifth.txt","rb") as fp:  # Pickling
        list_of_gene_objects = pickle.load(fp)

print('Adding exon information...')
gene_dict = Exon_Information.read_Ensembl_GRCh38_gtf_file_generate_nested_dict('/Users/jacob/Desktop/Isoform Mapper Webtool/Homo_sapiens.GRCh38_protein_coding.gtf')
genes_dict_median = Exon_Information.pick_exon_length_median_from_nested_dict(gene_dict)
Exon_Information.add_exon_median_to_gene_objects(list_of_gene_objects,genes_dict_median)
#Exon_Information.add_exon_objects_to_protein_objects(list_of_gene_objects,gene_dict)

print('Pickling list of gene objects and saving file...')
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_sixth.txt", "wb") as fp:  # Pickling
    pickle.dump(list_of_gene_objects, fp)

print('Validating and Correcting Library...') #missing variables and function input
print('Statistics before Clean-up:')
Validate_library.check_if_gene_name_and_prot_seq_are_switched(list_of_gene_objects)
gene_duplicates_dict = Validate_library.check_if_there_are_AA_seq_duplicates(list_of_gene_objects)[0]
print('Correcting library...')
list_of_gene_objects = Validate_library.fuse_attributes_of_duplicated_AA_seq_within_gene_object(list_of_gene_objects,gene_duplicates_dict)
list_of_gene_objects = Validate_library.delete_genes_and_protein_isoforms_with_no_or_one_AA_seq(list_of_gene_objects)
print('Statistics after Clean-up:')
Validate_library.check_if_gene_name_and_prot_seq_are_switched(list_of_gene_objects)
gene_duplicates_dict = Validate_library.check_if_there_are_AA_seq_duplicates(list_of_gene_objects)[0]
print('Pickling list of gene objects and saving file...')
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_25_march_final.txt", "wb") as fp:  # Pickling
    pickle.dump(list_of_gene_objects, fp)
print('Library Generation successfully executed.')
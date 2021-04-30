from Gene import *
from Protein_isoform import *
from Streamlit_community import *
from Input_flow import *
from Streamlit_Pop_ups import *
from Alignment import *
from Visualise_Alignment import *
from User_Input_Preparation import *
from Input_flow import *
from Table_Generation import *
from PIL import Image
from Statistics import *

#load library
list_of_gene_objects = Input_flow.import_data_from_github('list_of_gene_objects_19th_april.txt.gz')

#standard function parameters
match = 1
mismatch =-2
open_gap_penalty = -1.75
gap_extension_penalty = 0
minimal_exon_length = 5

#default include all columns
chosen_columns = ['Gene name', 'Ensembl Gene ID (ENSG)', 'Ensembl Transcript ID (ENST)',
                      'Ensembl Protein ID (ENSP)', 'Transcript name',
                      'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)', 'Refseq Protein ID (NP)',
                      'UCSC Stable ID (uc)',
                      'Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID', 'Uniparc ID',
                      'Ensembl Gene ID version (ENSG.Number)', 'Ensembl Transcript ID version (ENST.Number)',
                      'Ensembl Protein ID version (ENSP.Number)', 'Refseq Transcript ID version (NM.Number)',
                      'Refseq Transcript ID version (NP.Number)',
                      'HGNC ID (HGNC:Number)']
big_nested_dict = {}
for index,gene in enumerate(list_of_gene_objects):
    #create nested dict for one gene
    nested_dict = {gene.ensembl_gene_symbol:{index:Input_flow.pick_index_of_canonical_sequence(list_of_gene_objects,index)}}
    big_nested_dict.update(nested_dict)

#print(big_nested_dict)

#create big dataframe
big_df = Table_Generation.create_table_for_dict_of_gene_objects(big_nested_dict, list_of_gene_objects, chosen_columns,match, mismatch, open_gap_penalty,gap_extension_penalty)
big_df.to_csv('/Users/jacob/Desktop/all_genes_mapped.tsv')
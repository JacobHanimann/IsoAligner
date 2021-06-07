from Gene import *
from Protein_isoform import *
from Streamlit_community import *
from Input_flow import *
from Streamlit_Pop_ups import *
from Alignment import *
from Visualise_Alignment import *
from User_Input_Preparation import *
from Input_flow import *
from Table_Generation_match_statistics import *
from PIL import Image
from Statistics import *

#choose if standard needleman-wunsch should be used
conventional = False

#load library
list_of_gene_objects = Input_flow.import_data_from_github('../list_of_gene_objects_1st_june.txt.gz')

#standard function parameters
match = 1
mismatch =-2
open_gap_penalty = -1
gap_extension_penalty = 0

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
correct_aa_collection = 0
false_aa_collection = 0
mismatch_aa_collection = 0
for index,gene in enumerate(list_of_gene_objects):
    #create nested dict for one gene
    if len(gene.protein_sequence_isoform_collection)==1:
        continue
    if gene.minimal_exon_length!=None:
        nested_dict = {gene.ensembl_gene_symbol:{index:Input_flow.pick_index_of_canonical_sequence(list_of_gene_objects,index)}}
        big_nested_dict.update(nested_dict)
    if index % 1000 ==0:
        correct_aa, false_aa, mismatch_aa = Table_Generation_match.create_table_for_dict_of_gene_objects(
            big_nested_dict, list_of_gene_objects, chosen_columns, match, mismatch, open_gap_penalty,
            gap_extension_penalty, conventional=conventional)
        print(correct_aa)
        print(false_aa)
        print(mismatch_aa)
        correct_aa_collection = correct_aa_collection + correct_aa
        false_aa_collection = false_aa_collection + false_aa
        mismatch_aa_collection = mismatch_aa_collection + mismatch_aa
        nested_dict.clear()
        big_nested_dict.clear()


print(correct_aa_collection)
print(false_aa_collection)
print('mismatch',mismatch_aa_collection)
print(false_aa_collection/(correct_aa_collection+false_aa_collection))
from Ensembl import *
from HGNC import *
from Biomart_tables import *
from Refseq import *
from Uniprot import *
from Validation_of_library import *
from minimal_exon_length import *

date = '27th_may'




print('Pickling list of gene objects and saving file...')
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_"+date+"_fifth.txt", "rb") as fp:  # Pickling
    list_of_gene_objects = pickle.load(fp)


print('Adding exon information...')
gene_dict = Exon_Information.read_Ensembl_GRCh38_gtf_file_generate_nested_dict('/Users/jacob/Desktop/Isoform Mapper Webtool/Homo_sapiens.GRCh38_protein_coding.gtf')
genes_dict_median = Exon_Information.pick_exon_length_minimal_from_nested_dict(gene_dict)
Exon_Information.add_exon_minimal_to_gene_objects(list_of_gene_objects, genes_dict_median)
#Exon_Information.add_exon_objects_to_protein_objects(list_of_gene_objects,gene_dict)

count = 0
for gene in list_of_gene_objects:
    if gene.minimal_exon_length != None:
        count += 1

print(len(list_of_gene_objects))
print('genes with minimal exon length:', count)

gene_dict_2 = Exon_Information.read_Ensembl_GRCh38_gtf_file_generate_nested_dict('/Users/jacob/Desktop/GRCh38.103_protein_coding_patch.gtf')
genes_dict_median_2 = Exon_Information.pick_exon_length_minimal_from_nested_dict(gene_dict_2)
Exon_Information.add_exon_minimal_to_gene_objects(list_of_gene_objects, genes_dict_median_2)

count = 0
for gene in list_of_gene_objects:
    if gene.minimal_exon_length!=None:
        count +=1

print(len(list_of_gene_objects))
print('genes with minimal exon length second:', count)
import requests
from IsoAligner_core.Input_flow import *

list_of_gene_objects = Input_flow.import_data_from_github('Human_Isoform_Library/list_of_gene_objects_25th_july.txt.gz')

import requests
for gene in list_of_gene_objects:

    r= requests.get("http://www.isoaligner.org/api/map?id1="+gene.ensembl_gene_symbol)

print(r.text)
print(r.json())
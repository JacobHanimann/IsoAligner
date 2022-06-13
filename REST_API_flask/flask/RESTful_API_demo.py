import requests

def map_amino_acid_position_of_two_isoforms(reference_isoform_id, alternative_isoform_id):
    """query corresponding amino acid position between to isoforms"""
    #import python requests for API calls and json loading
    import json
    import requests
    #build query
    query = 'https://www.isoaligner.org/api/map?'+'id1='+ reference_isoform_id +'&id2='+ alternative_isoform_id
    #send API request
    response = requests.get(query, auth=False)
    #check if requests was successful
    print("Response status",response.status_code)
    #decode and load json format
    list_of_mapped_amino_acids= json.loads(response.json())
    return list_of_mapped_amino_acids


#EXAMPLE REQUEST

list_of_mapped_amino_acids = map_amino_acid_position_of_two_isoforms("EGFR-207", "EGFR-201")
#print(list_of_mapped_amino_acids)

#check if position 394 are the same for the two isoforms
for mapping in list_of_mapped_amino_acids:
    if mapping["ReferencePos"]==394:
        pass
        #print("Corresponding amino acid in alternative isoform:",mapping["IsoformPos"])


#INTEGRATE REQUESTS IN ANY ANALYSIS PIPELINE

#create a loop to check every isoform pair per gene in your dictionary
def compare_isoform_lists(gene_dict):
    for isoform1, isoform2 in gene_dict.keys():
        list_of_mapped_amino_acids = map_amino_acid_position_of_two_isoforms(isoform1, isoform2)
        pass

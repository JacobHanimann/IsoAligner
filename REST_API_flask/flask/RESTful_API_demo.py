import requests

def map_amino_acid_position_of_two_isoforms(reference_isoform_id, alternative_isoform_id):
    """query corresponding amino acid position between to isoforms"""
    #import python requests for API calls and json loading
    import json
    import requests
    #build query
    query = 'https://www.isoaligner.org/api/map?'+'id1='+ reference_isoform_id +'&id2='+ alternative_isoform_id
    #send api request
    response = requests.get(query, auth=False)
    #check if requests was successful
    print(response.status_code)
    #decode and load json format
    list_of_mapped_amino_acids= json.loads(response.json())
    return list_of_mapped_amino_acids


list_of_mapped_amino_acids = map_amino_acid_position_of_two_isoforms("EGFR-207", "EGFR-201")

# access position
print(list_of_mapped_amino_acids[300]["ReferencePos"])
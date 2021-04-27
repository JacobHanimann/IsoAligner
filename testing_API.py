import requests

BASE = "http://127.0.0.1:5000/"

data = {
 "sequence1": "NDFKLNDFDNFKLNSDLKFNASLKFNDSLJFNDSKF",
 "sequence2": "KDNFKLNFSDDNFKLNSDLKFNASLKFKNDFLKNS",
  "match": 2,
  "mismatch":4,
  "open_gap":3
}

response3 = requests.get(BASE+'map/positions?id1=BRAF-213&id2=BRAF-201&table_ids=[hgnc]&aa_position=49',data)

print(response3)
#print(response3.text)
print(response3.json())

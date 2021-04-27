import requests

BASE = "http://127.0.0.1:5000/"

data = {
 "sequence1": "NDFKLNDFDNFKLNSDLKFNASLKFNDSLJFNDSKF",
 "sequence2": "KDNFKLNFSDDNFKLNSDLKFNASLKFKNDFLKNS",
  "match": 2,
  "mismatch":4,
  "open_gap":3
}

response3 = requests.get(BASE+'map/positions?id1=KRAS-202&id2=KRAS-203&table_ids=[hgnc]',data)

print(response3)
#print(response3.text)
print(response3.json())

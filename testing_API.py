import requests

BASE = "http://127.0.0.1:5000/"

data = {
 "sequence1": "NDFKLNDFDNFKLNSDLKFNASLKFNDSLJFNDSKF",
  "match": 2,
  "mismatch":4,
  "open_gap":3,
}

for i in range (0,1):
 response3 = requests.get(BASE+'map?id1=EGFR',data)
 print(response3)
 print(response3.text)
 print(response3.json())
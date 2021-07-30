import requests

BASE = "http://jacobha.eu.pythonanywhere.com/"

data = {
 "seq1": "NDFKLNDFDNFKLNSDLKFNASLKFNDSLJFNDSKF",
"seq2": "NDFKLNDFDNFDFDKDKLNSDLKFNASLKFNDSLJFNDSKF",
  "match": 2,
  "mismatch":4,
  "open_gap":3,
}

for i in range (0,1):
 response3 = requests.get(BASE+'align?',data)
 print(response3)
 print(response3.text)
 print(response3.json())
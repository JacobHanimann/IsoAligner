import requests

BASE = "https://share.streamlit.io/jacobhanimann/isoaligner/main/Webinterface.py/"

data = {
 "sequence1": "NDFKLNDFDNFKLNSDLKFNASLKFNDSLJFNDSKF",
 "sequence2": "KDNFKLNFSDDNFKLNSDLKFNASLKFKNDFLKNS",
  "match": 2,
  "mismatch":4,
  "open_gap":3
}

response3 = requests.get(BASE+'map?id1=KRAS',data)

print(response3)
print(response3.text)
print(response3.json())
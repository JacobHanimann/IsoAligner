import requests

BASE = "http://127.0.0.1:5000/"

response = requests.get(BASE+'Align/ENSG58945')

print(response.json())
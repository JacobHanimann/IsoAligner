import requests


r= requests.get("http://www.isoaligner.org/api/map?id1=EGFR")

print(r.text)
print(r.json())
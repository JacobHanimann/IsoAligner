import requests


r= requests.get("https://isoaligner.org/api/map?")

print(r.text)
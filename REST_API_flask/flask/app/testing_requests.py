import requests


r= requests.get("https://isoaligner.org")

print(r.text)
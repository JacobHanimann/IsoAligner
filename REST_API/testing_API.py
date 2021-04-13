import requests

BASE = "http://127.0.0.1:5000/"

response = requests.put(BASE+'Align/DBFJDLJFNLDFBJFLDBDFNLK/DFBLJDFNLJFNLDFBJDBFLJBFDN')

print(response.json())

response2 = requests.get(BASE+'Convert/THisisID')

print(response2)
print(response2.json())




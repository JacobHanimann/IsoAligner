from flask import Flask

app = Flask(__name__)

@app.route("/")
def index():
    return "Hello from Isoaligner Flask"

if __name__ == "__main__":
    app.run(debug=True)

#from REST_API import *c
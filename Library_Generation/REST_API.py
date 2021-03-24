from flask import Flask
from flask_restful import Api, Resource


app = Flask(__name__)
api = Api(app)

class Helloworld(Resource):
    def get(self):
        return {'data':'Its me, mario'}

api.add_resource(Helloworld,'/HelloWorld')

if __name__ == "__main__":
    app.run(debug=True)
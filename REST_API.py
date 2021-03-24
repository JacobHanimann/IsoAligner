from flask import Flask
from flask_restful import Api, Resource


app = Flask(__name__)
api = Api(app)

class IsoAligner(Resource):
    def get(self,isoform_ID):
        return {'data':'IsoAlginer','isoform_ID_sent': isoform_ID}

api.add_resource(IsoAligner,'/Align/<string:isoform_ID>')

if __name__ == "__main__":
    app.run(debug=True)
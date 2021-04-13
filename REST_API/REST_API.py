from flask import Flask
from flask_restful import Api, Resource, reqparse, abort, fields, marshal_with
from flask_caching import Cache
from Visualise_Alignment import *
from Alignment import *
import gzip
import pickle

#Initialising Flask API and Cache
app = Flask(__name__)
api = Api(app)
cache = Cache()
app.config['CACHE_TYPE'] = 'simple'
cache.init_app(app)


data_args = reqparse.RequestParser()
data_args.add_argument("match", type=int, help="match score needed")
data_args.add_argument("mismatch", type=int, help="mismatch penalty needed",required=False)


@cache.cached(timeout=200,key_prefix='importing_library')
def import_data_from_github(file):
    '''import reference file (database), a pickle file generated in the database_old.py file'''
    with gzip.open(file, "rb") as fp:  # Pickling
        list_of_gene_objects = pickle.load(fp)
    return list_of_gene_objects

def align_sequences(input1,input2):
    needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(input1, input2, 1, -2,-1.75, 0,5)
    isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
    percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check, input1, input2)
    alignment_string =Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta, alignment_isoform_fasta,isoform_pattern_check, percentage_reference,percentage_isoform)
    return alignment_string


class IsoAligner(Resource):
    def get(self,isoform_ID):
        return {'data':'IsoAlginer','isoform_ID_sent': isoform_ID}

    def put(self,sequence1,sequence2):
        alignment= align_sequences(sequence1,sequence2)
        #list_of_gene_objects = import_data_from_github('../list_of_gene_objects_11_march.txt.gz')
        #return list_of_gene_objects[1].ensembl_gene_symbol
        return alignment

class Isoform(Resource):
    def post(self,number):
        args = data_args.parse_args()
        return args

data = {'isoform_ID':'name', 'gap_open':3}

api.add_resource(IsoAligner,'/Align/<string:sequence1>/<string:sequence2>','/Convert/<string:isoform_ID>')
api.add_resource(Isoform,'/Post/<int:number>')

if __name__ == "__main__":
    app.run(debug=True)
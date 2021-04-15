from flask import Flask
from flask_restful import Api, Resource, reqparse, abort, fields, marshal_with
from flask_caching import Cache
from Visualise_Alignment import *
from Alignment import *
import gzip
import pickle
from Gene import *
from Protein_isoform import *
import sys
# insert at position 1 in the path, as 0 is the path of this file.
#sys.path.insert(1, '../')

#Initialising Flask API and Cache
app = Flask(__name__)
api = Api(app)
cache = Cache()
app.config['CACHE_TYPE'] = 'simple'
cache.init_app(app)


# Arguments in the body of the requests
#mapping table
map_args = reqparse.RequestParser()
map_args.add_argument("match", type=int, help="set to default: 1", required=False)
map_args.add_argument("mismatch", type=int, help="set to default: -2", required=False)
map_args.add_argument("open_gap_penalty", type=float, help="set to default: -1.75", required=False)
map_args.add_argument("gap_extension_penalty", type=int, help="set to default: 0", required=False)
map_args.add_argument("minimal_exon_length", type=int, help="set to default: 5", required=False)
#raw alignment
align_args = reqparse.RequestParser()
align_args.add_argument("sequence1", type=str, help="reference raw amino acid required", required=True)
align_args.add_argument("sequence2", type=str, help="second raw amino acid required", required=True)
align_args.add_argument("match", type=int, help="set to default: 1", required=False)
align_args.add_argument("mismatch", type=int, help="set to default: -2", required=False)
align_args.add_argument("open_gap_penalty", type=float, help="set to default: -1.75", required=False)
align_args.add_argument("gap_extension_penalty", type=int, help="set to default: 0", required=False)
align_args.add_argument("minimal_exon_length", type=int, help="set to default: 5", required=False)


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


class Mapping_Table(Resource):
    def post(self, reference_ID, alternative="optional", aa_position='optional'):
        args = map_args.parse_args()
        list_of_gene_objects = import_data_from_github('list_of_gene_objects_11_march.txt.gz')
        if alternative!= 'optional':
            if aa_position!='optional':
                return list_of_gene_objects[1].ensembl_gene_symbol
            return alternative
        else:
            return reference_ID


class Raw_alignment(Resource):
    def post(self,option='mapping_table'):
        args = align_args.parse_args()
        return 'it worked'


api.add_resource(Mapping_Table,'/map/<string:reference_ID>','/map/<string:reference_ID>/<string:alternative>','/map/<string:reference_ID>/<string:alternative>/position/<int:aa_position>')
api.add_resource(Raw_alignment, '/align', '/align/visualise')

if __name__ == "__main__":
    app.run(debug=True)
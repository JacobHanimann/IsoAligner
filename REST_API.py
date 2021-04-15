from flask import Flask
from flask_restful import Api, Resource, reqparse, abort, fields, marshal_with
from flask_caching import Cache
from Visualise_Alignment import *
from Alignment import *
import pickle
from Gene import *
from Protein_isoform import *
import sys
# insert at position 1 in the path, as 0 is the path of this file.
#sys.path.insert(1, '../')
import sys
from Gene import *
from Protein_isoform import *
from Streamlit_community import *
from Input_flow import *
from Streamlit_Pop_ups import *
from Alignment import *
from Visualise_Alignment import *
from User_Input_Preparation import *
from Input_flow import *
from Table_Generation import *
from PIL import Image
from Statistics import *
from API_data_processing import *


#Initialising Flask API and Cache
app = Flask(__name__)
api = Api(app)
cache = Cache()
app.config['CACHE_TYPE'] = 'simple'
cache.init_app(app)

#importing isoform library
@cache.cached(timeout=300,key_prefix='importing_library') #makes no difference if function is cached or not
def import_data_from_github(file):
    '''import reference file (database), a pickle file generated in the database_old.py file'''
    with gzip.open(file, "rb") as fp:  # Pickling
        list_of_gene_objects = pickle.load(fp)
    return list_of_gene_objects

list_of_gene_objects = import_data_from_github('list_of_gene_objects_25th_march.txt.gz')

#standard parameters if no body is sent with the request
match = 1
mismatch =-2
open_gap_penalty = -1.75
gap_extension_penalty = 0
exon_length_AA= 5

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


def align_sequences(input1,input2):
    needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(input1, input2, 1, -2,-1.75, 0,5)
    isoform_pattern_check, alignment_reference_fasta, alignment_isoform_fasta = needleman_mapped[4:7]
    percentage_reference, percentage_isoform = Visualise_Alignment.calculate_percentage_of_mapped_positions(isoform_pattern_check, input1, input2)
    alignment_string =Visualise_Alignment.visualise_alignment_dynamically(alignment_reference_fasta, alignment_isoform_fasta,isoform_pattern_check, percentage_reference,percentage_isoform)
    return alignment_string


class Mapping_Table(Resource):
    def post(self, reference_ID, alternative="optional", aa_position='optional'):
        args = map_args.parse_args()
        if alternative!= 'optional':
            if aa_position!='optional':
                return list_of_gene_objects[200].ensembl_gene_symbol
        else: #just one reference ID given
            dict_of_IDs = Input_preparation.identify_IDs_from_user_text_input(reference_ID)
            gene_index = Input_flow.search_through_database_with_known_ID_Type(list_of_gene_objects, dict_of_IDs)
            verified_gene_index = Input_flow.remove_dict_elements_with_no_gene_object_match(gene_index)
            if not verified_gene_index:
                return 'ID not found'
            nested_dict = Input_flow.generate_nested_dictionary_with_index_of_canonical_protein_object(dict_of_IDs, verified_gene_index,list_of_gene_objects)
            #generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,
            #                                                                    list_of_gene_objects, index_gene_object,
            #                                                                    chosen_columns, match, mismatch,
            #                                                                    open_gap_penalty, gap_extension_penalty,
            #                                                                    exon_length_AA)
            return nested_dict


class Raw_alignment(Resource):
    def post(self,option='mapping_table'):
        args = align_args.parse_args()
        if option=='mapping_table':
            needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(args['sequence1'], args['sequence2'], match, mismatch,open_gap_penalty, gap_extension_penalty,exon_length_AA)
            generated_table = Table_Generation.create_pandas_dataframe_raw_aa_sequence(needleman_mapped)
            table_json = generated_table.to_json(orient='records')
            return table_json
        elif option=='visualise':
            alignment_string = Data_processing.align_sequences(args['sequence1'], args['sequence2'])
            return alignment_string


api.add_resource(Mapping_Table,'/map/<string:reference_ID>','/map/<string:reference_ID>/<string:alternative>','/map/<string:reference_ID>/<string:alternative>/position/<int:aa_position>')
api.add_resource(Raw_alignment, '/align', '/align/<string:option>')

if __name__ == "__main__":
    app.run(debug=True)
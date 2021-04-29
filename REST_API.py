from flask import Flask
from flask_restful import Api, Resource, reqparse, abort, fields, marshal_with
from flask_caching import Cache
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

list_of_gene_objects = import_data_from_github('list_of_gene_objects_19th_april.txt.gz')

#standard parameters if no body is sent with the request
match = 1
mismatch =-2
open_gap_penalty = -1.75
gap_extension_penalty = 0
exon_length_AA= 5
#standard ID's included in mapping table


# Arguments in the body of the requests
#mapping table
map_args = reqparse.RequestParser()
map_args.add_argument('id1',type=str, help="reference ID is required", required=True)
map_args.add_argument('id2',type=str, help="alternative ID is missing", required=False)
map_args.add_argument('aa_position',type=str, help="choose single AA position", required=False)
map_args.add_argument("match", type=int, help="set to default: 1", required=False)
map_args.add_argument("mismatch", type=int, help="set to default: -2", required=False)
map_args.add_argument("open_gap_penalty", type=float, help="set to default: -1.75", required=False)
map_args.add_argument("gap_extension_penalty", type=int, help="set to default: 0", required=False)
map_args.add_argument("minimal_exon_length", type=int, help="set to default: 5", required=False)
map_args.add_argument('table_ids', type=str, help='choose which IDs should be included in the mapping table',required=False)
#raw alignment
align_args = reqparse.RequestParser()
align_args.add_argument('view', type=bool, default=False, required=False, help='view alignment')
align_args.add_argument("sequence1", type=str, help="reference raw amino acid required", required=True)
align_args.add_argument("sequence2", type=str, help="second raw amino acid required", required=True)
align_args.add_argument("match", type=int, help="set to default: 1", required=False)
align_args.add_argument("mismatch", type=int, help="set to default: -2", required=False)
align_args.add_argument("open_gap_penalty", type=float, help="set to default: -1.75", required=False)
align_args.add_argument("gap_extension_penalty", type=int, help="set to default: 0", required=False)
align_args.add_argument("minimal_exon_length", type=int, help="set to default: 5", required=False)


#API methods
class Mapping_Table(Resource):
    def get(self):
        args = map_args.parse_args()
        if args['id2']!= None:
            nested_dict_reference,reference_type_dict = Data_processing.search_and_generate_nested_dict(args['id1'],list_of_gene_objects)
            nested_dict_alternative,alternative_type_dict = Data_processing.search_and_generate_nested_dict(args['id2'], list_of_gene_objects)
            if nested_dict_reference and nested_dict_alternative:
                index_of_gene = list(list(nested_dict_reference.values())[0].keys())[0]
                id_type_reference = list(reference_type_dict.values())[0]
                id_type_alternative = list(alternative_type_dict.values())[0]
                index_reference_transcript = list(list(nested_dict_reference.values())[0].values())[0]
                index_alternative_transcript = list(list(nested_dict_alternative.values())[0].values())[0]
                chosen_columns = Data_processing.choose_mapping_table_columns(args['table_ids'],id_type_reference,id_type_alternative,two_IDs=True)
                mapping_table = Data_processing.create_mapping_table_of_two_IDs(list_of_gene_objects, index_of_gene,
                                                                                index_reference_transcript,
                                                                                index_alternative_transcript,
                                                                                chosen_columns, match, mismatch,
                                                                                open_gap_penalty, gap_extension_penalty,
                                                                                exon_length_AA)
                if args['aa_position']==None:
                    table_json = mapping_table.to_json(orient='records')
                    return table_json
                else:
                    AA_new_position = Data_processing.extract_specific_position_mapping_table(mapping_table,args['aa_position'])
                    return AA_new_position
            if nested_dict_reference:
                return 'reference ID: '+args['id1']+' found, alternative ID: '+args['id2']+' not found', 400
            if nested_dict_alternative:
                return 'alternative ID: '+args['id2']+' found, reference ID: '+args['id1']+' not found', 400
            else:
                return 'IDs not found'

        else: #just one reference ID given
            nested_dict, dict = Data_processing.search_and_generate_nested_dict(args['id1'],list_of_gene_objects)
            if not nested_dict:
                return 'reference ID: '+args['id1']+' not found'
            index_gene_object = list(list(nested_dict.values())[0].keys())[0]
            index_of_reference_transcript = list(list(nested_dict.values())[0].values())[0]
            chosen_columns = Data_processing.choose_mapping_table_columns(table_ids=args['table_ids'])
            generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                list_of_gene_objects, index_gene_object,
                                                                                chosen_columns, match, mismatch,
                                                                                open_gap_penalty, gap_extension_penalty)
            table_json = generated_table.to_json(orient='records')
            return table_json


class Raw_alignment(Resource):
    def get(self):
        args = align_args.parse_args()
        if args['view']==False:
            needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(args['sequence1'], args['sequence2'], match, mismatch,open_gap_penalty, gap_extension_penalty,exon_length_AA)
            generated_table = Table_Generation.create_pandas_dataframe_raw_aa_sequence(needleman_mapped)
            table_json = generated_table.to_json(orient='records')
            return table_json
        elif args['view']==True:
            alignment_string = Data_processing.align_sequences(args['sequence1'], args['sequence2'])
            return alignment_string


#adding method to server
api.add_resource(Mapping_Table,'/map','/map/positions')
api.add_resource(Raw_alignment, '/align')

if __name__ == "__main__":
    app.run(debug=False)
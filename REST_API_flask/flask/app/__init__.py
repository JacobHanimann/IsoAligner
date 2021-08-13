from flask import Flask, render_template
from flask_restful import Api, Resource, reqparse
from cachetools import cached, TTLCache
from API_data_processing import *


#Initialising Flask API and Cache
app = Flask(__name__)
api = Api(app)
#cache = Cache(config={'CACHE_TYPE': 'SimpleCache'})
cache = TTLCache(maxsize=100, ttl=100)
app.config['CACHE_TYPE'] = 'simple'

@app.route("/")
def index():
    return render_template('iframe_streamlit.html')
#@cache.cached(timeout=300,key_prefix='importing_library') #makes no difference if function is cached or not
@app.route("/api")
def api_page():
    return "Go to REST API & Downloads on the main webpage to find out how to send requests"

@cached(cache)
def import_data_from_github(file):
    '''import reference file (database), a pickle file generated in the database_old.py file'''
    with gzip.open(file, "rb") as fp:  # Pickling
        list_of_gene_objects = pickle.load(fp)
    print('library loaded new')
    return list_of_gene_objects

#cache.set("list_of_gene_objects",import_data_from_github('list_of_gene_objects_4th_may.txt.gz'))
list_of_gene_objects = import_data_from_github('Human_Isoform_Library/list_of_gene_objects_25th_july.txt.gz')


# Arguments in the body of the requests
#mapping table
map_args = reqparse.RequestParser()
map_args.add_argument('id1',type=str, help="reference ID is required", required=True)
map_args.add_argument('id2',type=str, help="alternative ID is missing", required=False)
map_args.add_argument('pos',type=str, help="choose single AA position", required=False)
map_args.add_argument("match", type=int, help="set to default: 1, must be an integer", required=False)
map_args.add_argument("mismatch", type=int, help="set to default: -2, must be an integer", required=False)
map_args.add_argument("open_gap", type=float, help="set to default: -1", required=False)
map_args.add_argument("gap_ext", type=int, help="set to default: 0, must be an integer", required=False)
map_args.add_argument("min_ex_len", type=int, help="set to default: 12, must be an integer", required=False)
map_args.add_argument('df_ids', type=str, help='choose which IDs should be included in the mapping table',required=False)

#raw alignment
align_args = reqparse.RequestParser()
align_args.add_argument('view', type=bool, default=False, required=False, help='view alignment')
align_args.add_argument("seq1", type=str, help="reference raw amino acid required", required=True)
align_args.add_argument("seq2", type=str, help="second raw amino acid required", required=True)
align_args.add_argument("match", type=int, help="set to default: 1, must be an integer", required=False)
align_args.add_argument("mismatch", type=int, help="set to default: -2, must be an integer", required=False)
align_args.add_argument("open_gap", type=float, help="set to default: -1", required=False)
align_args.add_argument("gap_ext", type=int, help="set to default: 0, must be an integer", required=False)
align_args.add_argument("min_ex_len", type=int, help="set to default: 12, must be an integer", required=False)
align_args.add_argument("min_ex_len", type=int, help="set to default: 12", required=False)


#API methods
class Mapping_Table(Resource):
    def get(self):
        match = 1
        mismatch = -2
        open_gap_penalty = -1
        gap_extension_penalty = 0
        exon_length_AA = 12
        args = map_args.parse_args()
        if args['match'] != None:
            match = args['match']
        if args['mismatch'] != None:
            mismatch = args['mismatch']
        if args['open_gap'] != None:
            open_gap_penalty = args['open_gap']
        if args['gap_ext'] != None:
            gap_extension_penalty = args['gap_ext']
        if args['min_ex_len'] != None:
            exon_length_AA = args['min_ex_len']
        if args['id2']!= None:
            nested_dict_reference,reference_type_dict = Data_processing.search_and_generate_nested_dict(args['id1'],list_of_gene_objects)
            nested_dict_alternative,alternative_type_dict = Data_processing.search_and_generate_nested_dict(args['id2'], list_of_gene_objects)
            if nested_dict_reference and nested_dict_alternative:
                index_of_gene = list(list(nested_dict_reference.values())[0].keys())[0]
                id_type_reference = list(reference_type_dict.values())[0]
                id_type_alternative = list(alternative_type_dict.values())[0]
                index_reference_transcript = list(list(nested_dict_reference.values())[0].values())[0]
                index_alternative_transcript = list(list(nested_dict_alternative.values())[0].values())[0]
                chosen_columns = Data_processing.choose_mapping_table_columns(args['df_ids'],id_type_reference,id_type_alternative,two_IDs=True)
                mapping_table = Data_processing.create_mapping_table_of_two_IDs(list_of_gene_objects, index_of_gene,
                                                                                index_reference_transcript,
                                                                                index_alternative_transcript,
                                                                                chosen_columns, match, mismatch,
                                                                                open_gap_penalty, gap_extension_penalty,
                                                                                exon_length_AA)
                if args['pos']==None:
                    if type(mapping_table) != tuple:
                        table_json = mapping_table.to_json(orient='records')
                        return table_json
                    else:
                        return "No amino acid positions mapped. Adjust parameters to allow matches."
                else:
                    AA_new_position = Data_processing.extract_specific_position_mapping_table(mapping_table,args['pos'])
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
            chosen_columns = Data_processing.choose_mapping_table_columns(table_ids=args['df_ids'])

            generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                list_of_gene_objects, index_gene_object,
                                                                                chosen_columns, match, mismatch,
                                                                                open_gap_penalty, gap_extension_penalty)
            if type(generated_table) != tuple:
                table_json = generated_table.to_json(orient='records')
                return table_json
            else:
                return "No amino acid positions mapped. Adjust parameters to allow matches."


class Raw_alignment(Resource):
    def get(self):
        match = 1
        mismatch = -2
        open_gap_penalty = -1
        gap_extension_penalty = 0
        exon_length_AA = 12
        args = align_args.parse_args()
        if args['match'] != None:
            match = args['match']
        if args['mismatch'] != None:
            mismatch = args['mismatch']
        if args['open_gap'] != None:
            open_gap_penalty = args['open_gap']
        if args['gap_ext'] != None:
            gap_extension_penalty = args['gap_ext']
        if args['min_ex_len'] != None:
            exon_length_AA = args['min_ex_len']
        if args['view']==False:
            needleman_mapped = Alignment.map_AA_Needleman_Wunsch_with_exon_check(args['seq1'], args['seq2'], match, mismatch,open_gap_penalty, gap_extension_penalty,exon_length_AA)
            generated_table = Table_Generation.create_pandas_dataframe_raw_aa_sequence(needleman_mapped)
            if not generated_table.empty():
                table_json = generated_table.to_json(orient='records')
                return table_json
            else:
                return "No amino acid positions mapped. Adjust parameters to allow matches."
        elif args['view']==True:
            alignment_string = Data_processing.align_sequences(args['seq1'], args['seq2'])
            return alignment_string


#adding method
api.add_resource(Mapping_Table,'/api/map','/map/positions')
api.add_resource(Raw_alignment, '/api/align')
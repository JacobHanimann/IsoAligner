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

list_of_gene_objects = import_data_from_github('list_of_gene_objects_19th_april.txt.gz')

#standard parameters if no body is sent with the request
match = 1
mismatch =-2
open_gap_penalty = -1.75
gap_extension_penalty = 0
exon_length_AA= 5
#standard ID's included in mapping table
chosen_columns = ['Gene name', 'Ensembl Gene ID (ENSG)', 'Ensembl Transcript ID (ENST)', 'Ensembl Protein ID (ENSP)', 'Transcript name',
 'Refseq Gene ID (Number)', 'Refseq Transcript ID (NM)', 'Refseq Protein ID (NP)', 'UCSC Stable ID (uc)',
 'Uniprot Name ID', 'Uniprot Accession ID', 'Uniprot Isoform ID', 'Uniparc ID',
 'Ensembl Gene ID version (ENSG.Number)', 'Ensembl Transcript ID version (ENST.Number)',
 'Ensembl Protein ID version (ENSP.Number)', 'Refseq Transcript ID version (NM.Number)',
 'Refseq Transcript ID version (NP.Number)',
 'HGNC ID (HGNC:Number)']

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
    def post(self, reference_ID, alternative_ID="optional", aa_position='optional'):
        args = map_args.parse_args()
        if alternative_ID!= 'optional':
            if aa_position!='optional':
                return list_of_gene_objects[200].ensembl_gene_symbol
            nested_dict_reference = Data_processing.search_and_generate_nested_dict(reference_ID,list_of_gene_objects)
            nested_dict_alternative = Data_processing.search_and_generate_nested_dict(alternative_ID, list_of_gene_objects)
            if nested_dict_reference and nested_dict_alternative:
                index_of_gene = list(list(nested_dict_reference.values())[0].keys())[0]
                index_reference_transcript = list(list(nested_dict_reference.values())[0].values())[0]
                index_alternative_transcript = list(list(nested_dict_alternative.values())[0].values())[0]
                mapping_table = Data_processing.create_mapping_table_of_two_IDs(list_of_gene_objects,index_of_gene,index_reference_transcript,index_alternative_transcript,chosen_columns, match, mismatch,open_gap_penalty, gap_extension_penalty,exon_length_AA)
                table_json = mapping_table.to_json(orient='records')
                return table_json
            if nested_dict_reference:
                return 'reference ID found, alternative ID not found', 400
            if nested_dict_alternative:
                return 'alternative ID found, reference ID not found', 400
            else:
                return 'IDs not found'

        else: #just one reference ID given
            nested_dict = Data_processing.search_and_generate_nested_dict(reference_ID,list_of_gene_objects)
            if not nested_dict:
                return 'ID not found'
            index_gene_object = list(list(nested_dict.values())[0].keys())[0]
            index_of_reference_transcript = list(list(nested_dict.values())[0].values())[0]
            generated_table = Table_Generation.create_table_for_one_gene_object(index_of_reference_transcript,
                                                                                list_of_gene_objects, index_gene_object,
                                                                                chosen_columns, match, mismatch,
                                                                                open_gap_penalty, gap_extension_penalty,)
            table_json = generated_table.to_json(orient='records')
            return table_json


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


api.add_resource(Mapping_Table,'/map/<string:reference_ID>','/map/<string:reference_ID>/<string:alternative_ID>','/map/<string:reference_ID>/<string:alternative_ID>/position/<int:aa_position>')
api.add_resource(Raw_alignment, '/align', '/align/<string:option>')

if __name__ == "__main__":
    app.run(debug=True)
import streamlit as st
import sys
import SessionState
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
import requests


class Comparison():
    pass

    @staticmethod
    def extract_protein_sequences_in_fasta_from_gene(gene_object,Id_type="ENSP"):
        list_of_fasta_files = [">"+sequence.ENSP+"\n"+sequence.protein_sequence for sequence in gene_object.protein_sequence_isoform_collection if sequence.ENSP!=None]
        return list_of_fasta_files

    @staticmethod
    def post_fastas_to_clustal_MSA_API(fasta_string):
        baseUrl = u'https://www.ebi.ac.uk/Tools/services/rest/clustalo/'
        data = {'sequence': fasta_string,
                'outfmt': "fa",
                'stype': 'protein',
                'email': 'jacob.hanimann@gmail.com',
                'title': 'trialrun'}
        response3 = requests.post(baseUrl + 'run', data)
        return response3.text

    @staticmethod
    def get_MSA_alignment_clustal(fasta_string):
        baseUrl = u'https://www.ebi.ac.uk/Tools/services/rest/clustalo/'
        data = {'sequence': fasta_string,
                'outfmt': "fa",
                'stype': 'protein',
                'email': 'jacob.hanimann@gmail.com',
                'title': 'trialrun'}
        response3 = requests.get(baseUrl + 'run', data)
        return response3.text






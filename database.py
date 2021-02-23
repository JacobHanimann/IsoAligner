import pandas as pd
from collections.abc import Iterable
import pickle
from Gene import *
from Protein_isoform import *
from Alignment import *


def add_HCGN_information_to_gene_objects(file_of_gene_names,list_of_gene_objects):
    '''complement list of gene objects
    input: list of gene objects
    output: list of gene objects with added attribute values
    '''
    df = pd.read_csv(file_of_gene_names, sep='\t')
    print('total length: ',len(df))
    for index in range(0,len(df)):
        print(index)
        #extract data line by line
        HGNC = df.loc[index, 'HGNC']
        HGNC_gene_symbol = df.loc[index, 'approved_symbol']
        previous_symbols = df.loc[index, 'previous_symbols']
        refseq_gene_ID = df.loc[index, 'NCBI Gene ID']
        alias_symbols = df.loc[index, 'alias_symbols']
        uniprot_ID = df.loc[index, 'UniProt ID']  # check if it the same ID as in the Protein_isoform classes

        #transfrom data in correct format
        if type(previous_symbols) != float:  # None values are type float
            if "," in previous_symbols:
                previous_symbols = previous_symbols.split(', ')
            else:
                previous_symbols = [
                    previous_symbols]  # either way create a list because it facilitates later search functions
        if type(alias_symbols) != float:
            if "," in alias_symbols:
                alias_symbols = alias_symbols.split(', ')
            else:
                alias_symbols = [alias_symbols]

        found = False
        #search for a match
        for gene in list_of_gene_objects:
            if found:
                break
            if gene.ENSG == df.loc[index,'Ensembl gene ID']:
                  found = True
                  gene.HGNC = HGNC
                  gene.HGNC_gene_symbol = HGNC_gene_symbol
                  gene.previous_symbols = previous_symbols
                  gene.refseq_gene_ID = refseq_gene_ID
                  gene.alias_symbols = alias_symbols
                  gene.uniprot_ID = uniprot_ID #check if it the same ID as in the Protein_isoform classes

        if found == False:
            list_of_gene_objects.append(Gene(ENSG=df.loc[index,'Ensembl gene ID'],HGNC=HGNC, HGNC_gene_symbol = HGNC_gene_symbol, previous_symbols = previous_symbols, alias_symbols = alias_symbols, refseq_gene_ID=refseq_gene_ID))

    return list_of_gene_objects


def get_ensembl_fasta_sequences_and_IDs_and_create_gene_objects(file):
    '''extract fasta files one by one and add them to the gene objects'''
    with open(file, "r") as f:
        expenses_txt = f.readlines()
    # Put all the lines into a single string
    whole_txt = "".join(expenses_txt)
    splittext = re.split(">", whole_txt)
    fasta_count = 0
    matches = 0
    list_of_gene_objects =[]
    for fasta in splittext[1:]:
        fasta_count += 1
        found = False
        gene_name = get_bio_IDs_with_regex('gene_name', fasta)
        # create Protein_isoform object to add to the gene_object
        aa_sequence = Alignment.extract_only_AA_of_Fasta_file(fasta.split('\n', 1)[1])
        if aa_sequence==None: #sequence shorter than 7 AA long
            continue
        sequence_object = Protein_isoform(gene_name, aa_sequence,
                                          get_bio_IDs_with_regex('ensembl_ensg', fasta),
                                          get_bio_IDs_with_regex('ensembl_ensg_version', fasta),
                                          get_bio_IDs_with_regex('ensembl_enst', fasta),
                                          get_bio_IDs_with_regex('ensembl_enst_version', fasta),
                                          get_bio_IDs_with_regex('ensembl_ensp', fasta),
                                          get_bio_IDs_with_regex('ensembl_ensp_version', fasta),
                                          uniprot_accession=get_bio_IDs_with_regex('uniprot_accession', fasta),
                                          uniprot_uniparc=get_bio_IDs_with_regex('uniprot_uniparc', fasta))
        for gene in list_of_gene_objects:
            if found:
                break
            if gene.ENSG == sequence_object.ENSG:
                if type(gene.protein_sequence_isoform_collection)==list:
                    gene.protein_sequence_isoform_collection.append(sequence_object)
                else:
                    gene.protein_sequence_isoform_collection = [sequence_object]
                found=True

        if found ==False:
            list_of_gene_objects.append(Gene(sequence_object.ENSG,gene_name,protein_sequence_isoform_collection=[sequence_object]))

        print('Fasta files processed: ' + str(fasta_count) + '/' + str(len(splittext)))
    print('Fasta files matched: ' + str(matches))
    return list_of_gene_objects


def add_Uniprot_Isoform_refseqrna_transcript_name_ID_to_protein_attributes(file, list_of_gene_objects):
    '''
    add IDs to protein isoform object attributes to extend ID library
    :param file: from Biomart
    :return: updated list_of_gene_objects with more attributes
    '''
    df = pd.read_csv(file, sep='\t')
    print('total length:',len(df))
    for index in range(0,len(df)):
        print(index)
        found = False
        uniparc_ID = df.loc[index, 'UniParc ID']
        if type(uniparc_ID) ==float:
            continue
        for gene in list_of_gene_objects:
            if found:
                break
            if type(gene.protein_sequence_isoform_collection) == list:
                for sequence in gene.protein_sequence_isoform_collection:
                    if found:
                        break
                    if sequence.uniprot_uniparc == uniparc_ID:
                        found = True
                        if type(df.loc[index, 'RefSeq mRNA ID'])!=float:
                            sequence.refseq_NM = df.loc[index, 'RefSeq mRNA ID']
                        if type(df.loc[index, 'Transcript name'])!=float:
                            sequence.transcript_name = df.loc[index, 'Transcript name']
                        if type( df.loc[index, 'UniProtKB isoform ID'])!=float:
                            sequence.uniprot_isoform = df.loc[index, 'UniProtKB isoform ID']
                        break
            else:
                continue

def add_refseq_protein_IDs(file, list_of_gene_objects):
    '''add IDs from Biomart file'''
    df = pd.read_csv(file, sep='\t')
    print('total length',len(df))
    for index in range(0,len(df)):
        print(index)
        found = False
        uniparc_ID = df.loc[index, 'UniParc ID']
        if type(uniparc_ID)==float:
            continue
        for gene in list_of_gene_objects:
            if found:
                break
            if type(gene.protein_sequence_isoform_collection) == list:
                for sequence in gene.protein_sequence_isoform_collection:
                    if found:
                        break
                    if sequence.uniprot_uniparc == uniparc_ID:
                        found = True
                        if type(df.loc[index, 'RefSeq peptide ID'])!=float:
                            sequence.refseq_protein = df.loc[index, 'RefSeq peptide ID']
                            break
            else:
                continue

def add_refseq_fasta_sequences(file, list_of_gene_objects):
    '''complement library with fasta sequences from refseq
    idea: first check if IDs can be gene_found in the protein_isoform object
    if not: check if protein sequence is unique in the gene object
    if unique: add protein_isoform object
    if not unique: complement ID's'''

    #prepare file
    with open(file, "r") as f:
        expenses_txt = f.readlines()
        # Put all the lines into a single string
    whole_txt = "".join(expenses_txt)
    splittext = re.split("//\n", whole_txt)

    fasta_count = 0
    not_NP = 0
    HCGN_found = False
    NCBI_ID_found = False
    HCGN_count = 0
    NCBI_count = 0
    no_match =0
    sequences_added=0
    match_but_no_isoforms =0
    NP= False
    XP = False
    YP = False

    print('total length: ',len(splittext))
    for entry in splittext[0:-1]:
        fasta_count += 1
        print(fasta_count)

        #extract information out of entry
        #extract which type of ID is used
        try:
            Ids = re.findall("VERSION.*\n",entry)[0]
        except:
            print('did not find version')

        #extract exact IDs
        if "NP_" in Ids:
            NP_ID = get_bio_IDs_with_regex('refseq_prot',Ids)
            NP_version = get_bio_IDs_with_regex('refseq_prot_version',Ids)
            NM_ID_version = get_bio_IDs_with_regex('refseq_rna_version',re.findall("DBSOURCE.*\n",entry)[0])
            NP=True

        elif "XP_" in Ids:
            XP_ID = get_bio_IDs_with_regex('refseq_prot_predict', Ids)
            XP_version = get_bio_IDs_with_regex('refseq_prot_predict_version', Ids)
            XM_ID_version = get_bio_IDs_with_regex('refseq_rna_predict_version', re.findall("DBSOURCE.*\n", entry)[0])
            #print(entry)

        elif "YP_" in Ids:
            YP_ID = get_bio_IDs_with_regex('refseq_mitocho', Ids)
            YP_version = get_bio_IDs_with_regex('refseq_mitocho_version', Ids)
            NC_ID_version = get_bio_IDs_with_regex('refseq_chromosome_version', re.findall("DBSOURCE.*\n", entry)[0])
            print(entry)

        else:
            print('no known Version ID')
            print(Ids)
            not_NP += 1

        #try to extract HGNC_ID and NCBI_ID
        try:
            HGNC_ID = get_bio_IDs_with_regex('HGNC',re.findall('/db_xref="HGNC:HGNC:\d+',entry)[0])
            HCGN_found = True
            HCGN_count +=1
        except:
            pass
        try:
            NCBI_ID = re.findall('\d+',re.findall('/db_xref=\"GeneID:\d+\"',entry)[0])
            NCBI_ID_found = True
            NCBI_count += 1
        except:
            pass
        protein_sequence = extract_protein_sequence_from_refseq_entry(entry)

        #search for a match in gene list
        gene_found = False
        isoform_processed = False
        for gene in list_of_gene_objects:
            if isoform_processed:
                break
            if HCGN_found:
                if gene.HGNC==HGNC_ID:
                    gene_found = True
            if NCBI_ID_found:
                if gene.refseq_gene_ID == NCBI_ID:
                   gene_found = True
            if gene_found:
                if isoform_processed:
                    break
                #print(gene.HGNC, HGNC_ID)
                if NP:
                    if type(gene.protein_sequence_isoform_collection) == list:
                        for isoform in gene.protein_sequence_isoform_collection:
                            if isoform_processed:
                                break
                            #print(isoform.refseq_NP, NP_ID)
                            if isoform.refseq_NP == NP_ID:
                                if isoform.protein_sequence == protein_sequence:
                                    isoform.refseq_NP_version= NP_version
                                    isoform.refseq_NM_version=NM_ID_version
                                    print('ID versions added')
                                    isoform_processed = True
                                    break
                                else:
                                    print('same NP ID but not same sequence')
                            if isoform.protein_sequence == protein_sequence:
                                isoform.refseq_NP = NP_ID
                                isoform.refseq_NP_version = NP_version
                                isoform.refseq_NM_version = NM_ID_version
                                print('same sequence, added IDs know')

                        if not isoform_processed:
                                #print(gene.HGNC, HGNC_ID, NCBI_ID,gene.refseq_gene_ID,NP_ID)
                                gene.protein_sequence_isoform_collection.append(Protein_isoform(gene.ensembl_gene_symbol,protein_sequence,refseq_NM_version=NM_ID_version,refseq_NP=NP_ID, refseq_NP_version= NP_version))
                                isoform_processed = True
                                #print('added isoform')
                                sequences_added += 1
                elif XP:
                    if type(gene.protein_sequence_isoform_collection) == list:
                        gene.protein_sequence_isoform_collection.append(
                            Protein_isoform(gene.ensembl_gene_symbol, protein_sequence,
                                            refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,
                                            refseq_XP_version=XP_version))
                        isoform_processed = True
                        sequences_added += 1
                    else:
                        gene.protein_sequence_isoform_collection = [Protein_isoform(gene.ensembl_gene_symbol, protein_sequence,
                                            refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,
                                            refseq_XP_version=XP_version)]
                        isoform_processed = True
                        print('new collection XP')
                elif YP:
                    if type(gene.protein_sequence_isoform_collection) == list:
                        gene.protein_sequence_isoform_collection.append(
                        Protein_isoform(gene.ensembl_gene_symbol, protein_sequence,
                                        refseq_YP_version=YP_version, refseq_YP=YP_ID, refseq_NC_version=NC_ID_version))
                        isoform_processed = True
                        sequences_added += 1
                    else:
                        gene.protein_sequence_isoform_collection= [Protein_isoform(gene.ensembl_gene_symbol, protein_sequence,
                                        refseq_YP_version=YP_version, refseq_YP=YP_ID, refseq_NC_version=NC_ID_version)]
                        isoform_processed = True
                        print('new collection YP')

        #no match in list of gene objects
        #make a new gene object
        if not gene_found:
            if NP:
                if HCGN_found:
                    list_of_gene_objects.append(Gene(HGNC=HGNC_ID,refseq_gene_ID=NCBI_ID,
                                                     protein_sequence_isoform_collection=[Protein_isoform(protein_sequence,refseq_NM_version=NM_ID_version,refseq_NP=NP_ID,refseq_NP_version=NP_version)]))
                else:
                    list_of_gene_objects.append(Gene(refseq_gene_ID=NCBI_ID,
                                                     protein_sequence_isoform_collection=[Protein_isoform(protein_sequence,refseq_NM_version=NM_ID_version,refseq_NP=NP_ID,refseq_NP_version=NP_version)]))

            elif XP:
                if HCGN_found:
                    list_of_gene_objects.append(Gene(HGNC=HGNC_ID,refseq_gene_ID=NCBI_ID,
                                                     protein_sequence_isoform_collection=[Protein_isoform(protein_sequence,refseq_XM_version=XM_ID_version,refseq_XP=XP_ID,refseq_XP_version=XP_version)]))
                else:
                    list_of_gene_objects.append(Gene(refseq_gene_ID=NCBI_ID,
                                                     protein_sequence_isoform_collection=[
                                                         Protein_isoform(protein_sequence, refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,refseq_XP_version=XP_version)]))
            elif YP:
                if HCGN_found:
                    list_of_gene_objects.append(Gene(HGNC=HGNC_ID,refseq_gene_ID=NCBI_ID,
                                                     protein_sequence_isoform_collection=[Protein_isoform(protein_sequence,refseq_YP=YP_ID,refseq_NC_version=NC_ID_version,refseq_YP_version=YP_version)]))
                else:
                    list_of_gene_objects.append(Gene(refseq_gene_ID=NCBI_ID,
                                                     protein_sequence_isoform_collection=[
                                                         Protein_isoform(protein_sequence,refseq_YP=YP_ID,refseq_NC_version=NC_ID_version,refseq_YP_version=YP_version)]))


    print('total entries: ',len(splittext))
    print('not NP: ',not_NP)
    print('no matches: ',no_match)
    print('sequences added: ',sequences_added)
    print('match but no isoforms:',match_but_no_isoforms)
    print('HCGN count', HCGN_count)
    print('NCBI count', NCBI_count)


def extract_protein_sequence_from_refseq_entry(entry):
    '''extract protein sequence and format it correctly'''
    section = re.split("ORIGIN", entry)[1]
    without_newline = re.sub('\n', '', section)
    without_newline_and_whitespice = re.sub('\s', '', without_newline)
    without_numbers= re.sub("\d+",'',without_newline_and_whitespice)
    sequence_of_AA_acronym = '[A-Z]+'
    minimal_length_of_AA_seq = '[A-Z]{10}'
    raw_AA_seq_list = re.findall(sequence_of_AA_acronym + minimal_length_of_AA_seq, without_numbers.upper())
    if len(raw_AA_seq_list) >= 1:
        return raw_AA_seq_list[0]  # string
    else:
        return 'no AA sequence found'


def add_uniprot_fasta_files(file,list_of_objects):
    '''complement library with fasta sequences from uniprot'''


def get_bio_IDs_with_regex(ID_type, string):
    'generic functions to extract certain ID types from different databases'
    version = False
    #Ensembl
    if ID_type=='ensembl_ensg':
        pattern = 'ENSG\d{11}'
    elif ID_type == 'ensembl_ensg_version':
        pattern = 'ENSG\d+\.\d+'
        version = True
    elif ID_type == 'ensembl_enst':
         pattern = 'ENST\d{11}'
    elif ID_type == 'ensembl_enst_version':
         pattern = 'ENST\d+\.\d+'
         version = True
    elif ID_type == 'ensembl_ensp':
         pattern = 'ENSP\d{11}'
    elif ID_type == 'ensembl_ensp_version':
        pattern = 'ENSP\d+\.\d+'
        version = True

    #Refseq
    #known
    elif ID_type=='refseq_NM':
         pattern = 'NM_\d+'
    elif ID_type=='refseq_rna_version':
         pattern = 'NM_\d+\.\d+'
         version = True
    elif ID_type=='refseq_prot':
         pattern = 'NP_\d+'
    elif ID_type=='refseq_prot_version':
        version = True
        pattern = 'NP_\d+\.\d+'
    elif ID_type == 'refseq_prot_predict':
    #model refseq
        pattern = 'XP_\d+'
    elif ID_type == 'refseq_prot_predict_version':
        pattern = 'XP_\d+\.\d+'
        version = True
    elif ID_type == 'refseq_rna_predict':
        pattern = 'XM_\d+'
    elif ID_type == 'refseq_rna_predict_version':
        pattern = 'XM_\d+\.\d+'
        version = True
    #mitchondrium
    elif ID_type == 'refseq_mitocho':
        pattern = 'YP_\d+'
    elif ID_type == 'refseq_mitocho_version':
        version = True
        pattern = 'YP_\d+\.\d+'
    #refseq chromosome
    elif ID_type == 'refseq_chromosome':
        pattern = 'NC_\d+\.\d+'
    elif ID_type == 'refseq_chromosome_version':
        version = True
        pattern = 'NC_\d+\.\d+'

    #Uniprot IDs
    elif ID_type == 'uniprot_accession':
         pattern = '\|[OPQ][0-9][0-9A-Z]{3}[0-9]\||\|[A-NR-Z][0-9][A-Z][A-Z,0-9]{2}[0-9]\||\|[A-N,R-Z][0-9][A-Z][A-Z,0-9]{2}[0-9][A-Z][A-Z,0-9]{2}[0-9]\|'
    elif ID_type == 'uniprot_uniparc':
        pattern = 'UPI[0-9A-F]+'

    if ID_type == "gene_name":
        pattern = "\|[^\|\n]+\n"
        match_list = re.findall(pattern,string)
        if not match_list:  # if list is empty
            return 'not found'
        else:
            #remove \n with regex
            return match_list[0][1:-1] #remove \n

    #HGNC
    elif ID_type == 'HGNC':
        pattern = "HGNC:\d+"


#execute regular expression
    match_list = re.findall(pattern,string)
    if not match_list: #if list is empty
        return 'not found'
    elif len(match_list)==1:
        if ID_type == "uniprot_accession":
            return match_list[0][1:-1]
        else:
            return match_list[0]
    else:
        if version == False:
            return match_list[0]
        else:
            return match_list



def save_all_data_in_pickle_style():
    'write general function to save files in the pickle format'


def save_results_to_tsv_file(dictionary):
    'to be pre-computed values'



#Execution

##create list of gene objects
#print('generating gene list')
#list_of_gene_objects = get_ensembl_fasta_sequences_and_IDs_and_create_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/ensembl_fasta_IDs_gene_name.txt')
#
##add HCGN adn NCBI gene information
#print('adding HCGN information')
#add_HCGN_information_to_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/HGNC_protein_coding_ensembl.txt',list_of_gene_objects)
#
##add ID's to protein_isoform class
#print('add uniprot IDs')
#add_Uniprot_Isoform_refseqrna_transcript_name_ID_to_protein_attributes('/Users/jacob/Desktop/Isoform Mapper Webtool/NM_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
#print('add refseq Ids')
#add_refseq_protein_IDs('/Users/jacob/Desktop/Isoform Mapper Webtool/NP_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
#
##save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_22_feb.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)

with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_22_feb.txt", "rb") as fp:  # Pickling
        list_of_gene_objects = pickle.load(fp)

for gene in list_of_gene_objects:
    if type(gene.protein_sequence_isoform_collection)==list:
        for isoform in gene.protein_sequence_isoform_collection:
            if isoform.refseq_NM!=None:
                print(isoform.refseq_NM)

add_refseq_protein_IDs('/Users/jacob/Desktop/Isoform Mapper Webtool/NP_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)


#save list of gene objects to import to the subsequent script
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_23_feb_second.txt", "wb") as fp:  # Pickling
    pickle.dump(list_of_gene_objects, fp)


#add refseq fasta files and IDs
print('add refseq fasta')
add_refseq_fasta_sequences('/Users/jacob/Desktop/Isoform Mapper Webtool/refseq_fasta_and_info/GCF_000001405.39_GRCh38.p13_protein.gpff',list_of_gene_objects)

#save list of gene objects to import to the subsequent script
with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_23_feb_third.txt", "wb") as fp:  # Pickling
    pickle.dump(list_of_gene_objects, fp)

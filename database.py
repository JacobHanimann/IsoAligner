import pandas as pd
from collections.abc import Iterable
import pickle
from Gene import *
from Protein_isoform import *
from Alignment import *
import streamlit as st


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
        sequence_object = Protein_isoform(aa_sequence, gene_name,
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
    not_any_type = 0
    HCGN_count = 0
    NCBI_count = 0
    no_match =0
    sequences_added=0
    match_but_no_isoforms =0
    NP_ID_not_same_sequence = 0
    NP_ID_version_added = 0
    same_sequence_added_IDs = 0
    XP_added = 0
    YP_added = 0


    print('total length: ',len(splittext))
    for entry in splittext[0:-1]:

        #counting info
        NP = False
        XP = False
        YP = False
        HCGN_found = False
        NCBI_ID_found = False
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
            XP = True
            #print(entry)

        elif "YP_" in Ids:
            YP_ID = get_bio_IDs_with_regex('refseq_mitocho', Ids)
            YP_version = get_bio_IDs_with_regex('refseq_mitocho_version', Ids)
            NC_ID_version = get_bio_IDs_with_regex('refseq_chromosome_version', re.findall("DBSOURCE.*\n", entry)[0])
            YP = True
            #print(entry)

        else:
            print('no known Version ID')
            print(Ids)
            not_any_type += 1

        #try to extract HGNC_ID and NCBI_ID
        try:
            HGNC_ID = get_bio_IDs_with_regex('HGNC',re.findall('/db_xref="HGNC:HGNC:\d+',entry)[0])
            HCGN_found = True
            HCGN_count +=1
        except:
            pass
        try:
            NCBI_ID = int(re.findall('\d+',re.findall('/db_xref=\"GeneID:\d+\"',entry)[0])[0])
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
                                    NP_ID_version_added +=1
                                    isoform_processed = True
                                    break
                                else:
                                    print('same NP ID but not same sequence')
                                    NP_ID_not_same_sequence += 1
                            elif isoform.protein_sequence == protein_sequence:
                                isoform.refseq_NP = NP_ID
                                isoform.refseq_NP_version = NP_version
                                isoform.refseq_NM_version = NM_ID_version
                                same_sequence_added_IDs +=1
                                isoform_processed = True

                        if not isoform_processed:
                                #print(gene.HGNC, HGNC_ID, NCBI_ID,gene.refseq_gene_ID,NP_ID)
                                gene.protein_sequence_isoform_collection.append(Protein_isoform(protein_sequence,gene_name=gene.ensembl_gene_symbol,refseq_NM_version=NM_ID_version,refseq_NP=NP_ID, refseq_NP_version= NP_version))
                                isoform_processed = True
                                #print('added isoform')
                                sequences_added += 1
                elif XP:
                    XP_added += 1
                    isoform_processed = True
                    if type(gene.protein_sequence_isoform_collection) == list:
                        gene.protein_sequence_isoform_collection.append(
                            Protein_isoform(protein_sequence,gene_name=gene.ensembl_gene_symbol,
                                            refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,
                                            refseq_XP_version=XP_version))
                    else:
                        gene.protein_sequence_isoform_collection = [Protein_isoform(protein_sequence,gene_name=gene.ensembl_gene_symbol,
                                            refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,
                                            refseq_XP_version=XP_version)]
                        print('new collection XP')
                elif YP:
                    YP_added += 1
                    isoform_processed = True
                    if type(gene.protein_sequence_isoform_collection) == list:
                        gene.protein_sequence_isoform_collection.append(
                        Protein_isoform( protein_sequence,gene_name=gene.ensembl_gene_symbol,
                                        refseq_YP_version=YP_version, refseq_YP=YP_ID, refseq_NC_version=NC_ID_version))
                    else:
                        gene.protein_sequence_isoform_collection= [Protein_isoform(protein_sequence, gene_name=gene.ensembl_gene_symbol,
                                        refseq_YP_version=YP_version, refseq_YP=YP_ID, refseq_NC_version=NC_ID_version)]
                        print('new collection YP')

        #no match in list of gene objects
        #make a new gene object
        if not gene_found:
            isoform_processed=True
            no_match += 1
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
    print('not any type: ',not_any_type)
    print('no matches: ',no_match)
    print('NP sequences added: ',sequences_added)
    print('match but no isoforms:',match_but_no_isoforms)
    print('HCGN count', HCGN_count)
    print('NCBI count', NCBI_count)
    print('NP_ID_not_same_Seq',NP_ID_not_same_sequence)
    print('NP_Version added',NP_ID_version_added)
    print('same sequence found, added IDs',same_sequence_added_IDs)
    print('XP added',XP_added)
    print('YP added',YP_added)


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

    #prepare file
    with open(file, "r") as f:
        expenses_txt = f.readlines()
        # Put all the lines into a single string
    whole_txt = "".join(expenses_txt)
    splittext = re.split("\n>", whole_txt)

    #organisation
    no_gene_name = 0
    already_in_accession = 0
    accession_in_but_other_sequence = 0
    already_in_uniprot_isoform = 0
    uniprot_isoform_not_same_seq = 0
    new_isoform_for_gene = 0
    no_gene_match_found = 0
    fasta_count = 0
    protein_sequence_already_without_uniprot_ID = 0


    print(len(splittext))

    #iterate
    for fasta in splittext[1:len(splittext)]:

        #organisation
        gene_name_found = True
        fasta_count += 1
        print(fasta_count)
        uniprot_isoform = False

        #extracting information
        try:
            accession = re.split('\|',fasta)[1]
            if "-" in accession:
                uniprot_isoform= True
        except:
            print('no accession number')
            print(fasta)
        try:
            gene_name = re.findall("GN=[A-Z,0-9]+",fasta)[0][3:]
        except:
            no_gene_name +=1
            gene_name_found=False
        try:
            protein_sequence = Alignment.extract_only_AA_of_Fasta_file(re.split("\n",fasta,maxsplit=1)[1])
        except:
            print('no AA sequence found')

        #if fasta file contains a gene name
        if gene_name_found:

            #find gene object
            gene_identified = False
            for index, gene in enumerate(list_of_gene_objects):
                if gene_identified:
                    break
                list_of_gene_names = [gene.ensembl_gene_symbol] + [gene.HGNC_gene_symbol]
                if type(gene.alias_symbols) ==list:
                    list_of_gene_names = list_of_gene_names + gene.alias_symbols
                if type(gene.previous_symbols)==list:
                    list_of_gene_names = list_of_gene_names + gene.previous_symbols
                if gene_name in list_of_gene_names:
                    gene_identified= True
                    gene_index = index
                    break

            # if gene name was found in list of gene objects
            if gene_identified:
                if type(gene.protein_sequence_isoform_collection)==list:
                    found = False
                    for isoform in gene.protein_sequence_isoform_collection:
                        if found:
                            break
                        if isoform.uniprot_accession == accession:
                            if isoform.protein_sequence == protein_sequence:
                                already_in_accession += 1
                                isoform.uniprot_isoform = accession+'-1'
                                found = True
                            else:
                                accession_in_but_other_sequence += 1
                                isoform.uniprot_accession = None  # delete (false) attribute of isoform
                                gene.protein_sequence_isoform_collection.append(Protein_isoform(protein_sequence, uniprot_accession=accession, uniprot_isoform=accession+"-1", gene_name=gene_name)) #add isoform to collection
                                found = True

                        elif isoform.uniprot_isoform == accession:
                            if isoform.protein_sequence == protein_sequence:
                                already_in_uniprot_isoform += 1
                                found = True
                            else:
                                uniprot_isoform_not_same_seq +=1
                                isoform.uniprot_isoform = None  # delete (false) attribute of isoform
                                gene.protein_sequence_isoform_collection.append(
                                    Protein_isoform(protein_sequence, uniprot_accession=get_bio_IDs_with_regex('uniprot_accession',accession),uniprot_isoform=accession, gene_name=gene_name))
                                found = True

                        elif isoform.protein_sequence == protein_sequence:
                             found = True
                             protein_sequence_already_without_uniprot_ID +=1
                             if uniprot_isoform:
                                isoform.uniprot_accession=get_bio_IDs_with_regex('uniprot_accession',accession)
                                isoform.uniprot_isoform = accession
                             else:
                                 isoform.uniprot_accession = accession
                                 isoform.uniprot_isoform = accession+'-1'

                    if not found:
                        new_isoform_for_gene += 1
                        if uniprot_isoform:
                            gene.protein_sequence_isoform_collection.append(Protein_isoform(protein_sequence,uniprot_accession=get_bio_IDs_with_regex('uniprot_accession', accession),
                                            uniprot_isoform=accession, gene_name=gene_name))
                        else:
                            gene.protein_sequence_isoform_collection.append(Protein_isoform(protein_sequence, uniprot_accession=accession,uniprot_isoform=accession+'-1', gene_name=gene_name))



            # gene name was not found in list of gene objects
            else:
                no_gene_match_found +=1
                pass


    print('fasta files with no gene names:', no_gene_name)
    print('already in accession',already_in_accession)
    print('accession but other sequence',accession_in_but_other_sequence)
    print('already in uniprot isoform',already_in_uniprot_isoform)
    print('uniprot_isoform_not_same_seq',uniprot_isoform_not_same_seq)
    print('new_isoform_for_gene',new_isoform_for_gene)
    print('no_gene_match_found',no_gene_match_found)
    print('protein sequence match but no uniprot_ID associated: ', protein_sequence_already_without_uniprot_ID)


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


def check_if_there_are_AA_seq_duplicates(list_of_gene_objects):
    '''
    check out if there were IDs and Seq that are the same but escaped the match
    :param list_of_gene_objects:
    :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
    '''

    def duplicates(lst, item):
        return [i for i, x in enumerate(lst) if x == item]

    genes_without_AA_seq = 0
    duplicates_number = 0
    genes_without_duplicates = 0
    genes_with_more_than_one_duplicate =0
    redundant_sequences = 0
    duplicate_genes_dict = dict()
    for index,gene in enumerate(list_of_gene_objects):
        if type(gene.protein_sequence_isoform_collection)==list:
            List = [sequence.protein_sequence for sequence in gene.protein_sequence_isoform_collection]
            duplicates_dict = dict((x, duplicates(List, x)) for x in set(List) if List.count(x) > 1)
            if len(duplicates_dict) !=0:
                if list(duplicates_dict.keys())[0]!=None:
                    duplicate_genes_dict[index] = duplicates_dict
                    duplicates_number += 1
                    for sequence,objects in duplicates_dict.items():
                        redundant_sequences = redundant_sequences + len(objects)
                if len(duplicates_dict) >1:
                    genes_with_more_than_one_duplicate +=1
            else:
                genes_without_duplicates +=1
        else:
            genes_without_AA_seq += 1

    print('number of genes: ', len(list_of_gene_objects))
    print('genes with no AA seq: ', genes_without_AA_seq)
    print('number of genes with AA seq duplicates: ',duplicates_number)
    print('number of genes without AA seq duplicates: ',genes_without_duplicates)
    print('number of redundant AA sequences:', redundant_sequences)
    print('number of genes with more than one AA seq duplicate: ',genes_with_more_than_one_duplicate)
    return duplicate_genes_dict, duplicates_number,genes_without_duplicates,redundant_sequences,genes_with_more_than_one_duplicate


def check_if_there_are_exact_duplicates(list_of_gene_objects):
    '''
    check out if there were IDs and Seq that are the same but escaped the match
    :param list_of_gene_objects:
    :return: dictionary of gene indexes as key and dictionary of duplicates of the isoform collection
    '''

    def duplicates(lst, item):
        return [i for i, x in enumerate(lst) if x == item]

    genes_without_AA_seq = 0
    duplicates_number = 0
    genes_without_duplicates = 0
    genes_with_more_than_one_duplicate =0
    duplicate_genes_dict = dict()
    for index,gene in enumerate(list_of_gene_objects):
        if type(gene.protein_sequence_isoform_collection)==list:
            List = [tuple(list(sequence.__dict__.items())) for sequence in gene.protein_sequence_isoform_collection]
            duplicates_dict = dict((x, duplicates(List, x)) for x in set(List) if List.count(x) > 1)
            if len(duplicates_dict) !=0:
                if list(duplicates_dict.keys())[0]!=None:
                    duplicate_genes_dict[index] = duplicates_dict
                    duplicates_number += 1
                if len(duplicates_dict) >1:
                    genes_with_more_than_one_duplicate +=1
            else:
                genes_without_duplicates +=1
        else:
            genes_without_AA_seq += 1

    print('number of genes: ', len(list_of_gene_objects))
    print('genes with no AA seq: ', genes_without_AA_seq)
    print('number of genes with exact duplicates: ',duplicates_number)
    print('number of genes without exact duplicates: ',genes_without_duplicates)
    print('number of genes with more than one exact duplicate: ',genes_with_more_than_one_duplicate)
    return duplicate_genes_dict



def fuse_attributes_of_duplicated_AA_seq_within_gene_object(list_of_gene_objects,duplicate_genes_dict):
    '''
    function that fuses protein isoform objects if the attributes can complement each other to one big object, otherwise the duplicates will stay separated.
    :param list_of_gene_objects:
    :param duplicate_genes_dict:
    :return: updated list_of_gene_objects
    '''
    reduced_isoform_count = 0
    couldnotmatch = 0
    duplicates_in_total = 0
    for gene,duplicates_dict in duplicate_genes_dict.items():
        tobedeleted= []
        duplicates_in_total = duplicates_in_total + len(duplicates_dict)
        for duplicate_AA in duplicates_dict.items():
            new_object_attributes = Protein_isoform(duplicate_AA[0])
            isoform_dict = dict()
            list_of_attributes = [a for a in dir(new_object_attributes) if not a.startswith('__')]
            different_attributes = False
            for isoform in duplicate_AA[1]:
                if different_attributes:
                    break
                isoform = list_of_gene_objects[gene].protein_sequence_isoform_collection[isoform]
                for attribute in list_of_attributes:
                    if different_attributes:
                        break
                    if getattr(new_object_attributes,attribute)==None:
                        if getattr(isoform,attribute)!=None:
                            setattr(new_object_attributes,attribute,getattr(isoform,attribute))
                    else:
                        if getattr(isoform, attribute) != None:
                           if getattr(isoform, attribute) == getattr(new_object_attributes,attribute):
                               pass #attributes are the same, protein object can still be fused
                           else: #stop process, IDs differ from each other
                               different_attributes=True
                               couldnotmatch +=1
            if different_attributes== False:
                tobedeleted.extend(duplicate_AA[1])
                list_of_gene_objects[gene].protein_sequence_isoform_collection.append(new_object_attributes)
                reduced_isoform_count +=1
        if tobedeleted:
            for ele in sorted(tobedeleted, reverse=True):
                del list_of_gene_objects[gene].protein_sequence_isoform_collection[ele]

    print('duplicates in total:', duplicates_in_total)
    print('duplicates that could not be matched:',couldnotmatch)
    print('duplicates that could be matched:',reduced_isoform_count)
    return list_of_gene_objects


def check_if_gene_name_and_prot_seq_are_switched(list_of_gene_objects):
    '''somewhere in the database generation gene name and protein sequence attribute of a protein isoform object are being falsely switched'''
    false_assigned_gene_name_isoform = 0
    for gene in list_of_gene_objects:
        if type(gene.protein_sequence_isoform_collection) == list:
            for isoform in gene.protein_sequence_isoform_collection:
                if type(isoform.gene_name)==str:
                    if Alignment.extract_only_AA_of_Fasta_file(isoform.gene_name)!= None:
                        false_assigned_gene_name_isoform +=1
    print('number of falsely assigned AA seq to gene_name:',false_assigned_gene_name_isoform)


#Execution

##create list of gene objects
#print('generating gene list')
#list_of_gene_objects = get_ensembl_fasta_sequences_and_IDs_and_create_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/ensembl_fasta_IDs_gene_name.txt')
#
##save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_24_feb_first.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
##add HCGN adn NCBI gene information
#print('adding HCGN information')
#add_HCGN_information_to_gene_objects('/Users/jacob/Desktop/Isoform Mapper Webtool/HGNC_protein_coding_ensembl.txt',list_of_gene_objects)
#
##save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_24_feb_second.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
##add ID's to protein_isoform class
#print('add uniprot IDs')
#add_Uniprot_Isoform_refseqrna_transcript_name_ID_to_protein_attributes('/Users/jacob/Desktop/Isoform Mapper Webtool/NM_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
#
##save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_24_feb_third.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
#print('add refseq Ids')
#add_refseq_protein_IDs('/Users/jacob/Desktop/Isoform Mapper Webtool/NP_Uniprot_Isoform_uniparc.txt',list_of_gene_objects)
#
##save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_24_feb_fourth.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
##add refseq fasta files and IDs
#print('add refseq fasta')
#add_refseq_fasta_sequences('/Users/jacob/Desktop/Isoform Mapper Webtool/refseq_fasta_and_info/GCF_000001405.39_GRCh38.p13_protein.gpff',list_of_gene_objects)
#
##save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_24_feb_fifth.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)


#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_24_feb_fourth.txt", "rb") as fp:  # Pickling
#        list_of_gene_objects = pickle.load(fp)
#
#gene_duplicates_dict =check_if_there_are_exact_duplicates(list_of_gene_objects)
#
#fuse_duplicated_AA_seq_within_gene_object(list_of_gene_objects,gene_duplicates_dict)
#
#add_refseq_fasta_sequences('/Users/jacob/Desktop/Isoform Mapper Webtool/refseq_fasta_and_info/GCF_000001405.39_GRCh38.p13_protein.gpff',list_of_gene_objects)
#
###save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_9_march_first.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)
#
#gene_duplicates_dict =check_if_there_are_AA_seq_duplicates(list_of_gene_objects)
#
#fuse_duplicated_AA_seq_within_gene_object(list_of_gene_objects,gene_duplicates_dict)
#
#add_uniprot_fasta_files('/Users/jacob/Desktop/Isoform Mapper Webtool/uniprot_downloads/uniprot-proteome_UP000005640.fasta',list_of_gene_objects)
#
###save list of gene objects to import to the subsequent script
#with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_9_march_second.txt", "wb") as fp:  # Pickling
#    pickle.dump(list_of_gene_objects, fp)

#for gene in list_of_gene_objects:
#    if type(gene.protein_sequence_isoform_collection)==list:
#        for isoform in gene.protein_sequence_isoform_collection:
#            if isoform.refseq_NM!=None:
#                print(isoform.refseq_NM)

with open("/Users/jacob/Desktop/Isoform Mapper Webtool/list_of_gene_objects_with_fasta_9_march_second.txt","rb") as fp:  # Pickling
        list_of_gene_objects = pickle.load(fp)

check_if_there_are_exact_duplicates(list_of_gene_objects)

check_if_gene_name_and_prot_seq_are_switched(list_of_gene_objects)

gene_duplicates_dict =check_if_there_are_AA_seq_duplicates(list_of_gene_objects)[0]

list_of_gene_objects = fuse_attributes_of_duplicated_AA_seq_within_gene_object(list_of_gene_objects,gene_duplicates_dict)

print('UDPATED')

check_if_there_are_exact_duplicates(list_of_gene_objects)
check_if_gene_name_and_prot_seq_are_switched(list_of_gene_objects)
gene_duplicates_dict =check_if_there_are_AA_seq_duplicates(list_of_gene_objects)[0]
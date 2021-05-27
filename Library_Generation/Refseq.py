from Gene import *
from Protein_isoform import *
from Alignment import *
from Extractions_BioIDs import *


class Refseq():
    pass


    @staticmethod
    def add_refseq_fasta_sequences(file, list_of_gene_objects):
        '''complement library with fasta sequences from refseq
        idea: first check if IDs can be gene_found in the protein_isoform object
        if not: check if protein sequence is unique in the gene object
        if unique: add protein_isoform object
        if not unique: complement ID's'''

        # prepare file
        with open(file, "r") as f:
            expenses_txt = f.readlines()
            # Put all the lines into a single string
        whole_txt = "".join(expenses_txt)
        splittext = re.split("//\n", whole_txt)

        fasta_count = 0
        not_any_type = 0
        HCGN_count = 0
        NCBI_count = 0
        no_match = 0
        sequences_added = 0
        match_but_no_isoforms = 0
        NP_ID_not_same_sequence = 0
        NP_ID_version_added = 0
        same_sequence_added_IDs = 0
        XP_added = 0
        YP_added = 0

        print('total length: ', len(splittext))
        for entry in splittext[0:-1]:

            # counting info
            NP = False
            XP = False
            YP = False
            HCGN_found = False
            NCBI_ID_found = False
            fasta_count += 1
            if fasta_count % 1000 == 0:
                print(100 * round(fasta_count / len(splittext), 2), '%')

            # extract information out of entry
            # extract which type of ID is used
            try:
                Ids = re.findall("VERSION.*\n", entry)[0]
            except:
                print('did not find version')

            # extract exact IDs
            if "NP_" in Ids:
                NP_ID =Get_Bio_ID.get_bio_IDs_with_regex('refseq_prot', Ids)
                NP_version =Get_Bio_ID.get_bio_IDs_with_regex('refseq_prot_version', Ids)
                NM_ID_version =Get_Bio_ID.get_bio_IDs_with_regex('refseq_rna_version', re.findall("DBSOURCE.*\n", entry)[0])
                NP = True

            elif "XP_" in Ids:
                XP_ID =Get_Bio_ID.get_bio_IDs_with_regex('refseq_prot_predict', Ids)
                XP_version =Get_Bio_ID.get_bio_IDs_with_regex('refseq_prot_predict_version', Ids)
                XM_ID_version =Get_Bio_ID.get_bio_IDs_with_regex('refseq_rna_predict_version',
                                                       re.findall("DBSOURCE.*\n", entry)[0])
                XP = True
                # print(entry)

            elif "YP_" in Ids:
                YP_ID =Get_Bio_ID.get_bio_IDs_with_regex('refseq_mitocho', Ids)
                YP_version =Get_Bio_ID.get_bio_IDs_with_regex('refseq_mitocho_version', Ids)
                NC_ID_version =Get_Bio_ID.get_bio_IDs_with_regex('refseq_chromosome_version',
                                                       re.findall("DBSOURCE.*\n", entry)[0])
                YP = True
                # print(entry)

            else:
                print('no known Version ID')
                print(Ids)
                not_any_type += 1

            # try to extract HGNC_ID and NCBI_ID
            try:
                HGNC_ID =Get_Bio_ID.get_bio_IDs_with_regex('HGNC', re.findall('/db_xref="HGNC:HGNC:\d+', entry)[0])
                HCGN_found = True
                HCGN_count += 1
            except:
                pass
            try:
                NCBI_ID = str(re.findall('\d+', re.findall('/db_xref=\"GeneID:\d+\"', entry)[0])[0])
                NCBI_ID_found = True
                NCBI_count += 1
            except:
                pass
            protein_sequence = Refseq.extract_protein_sequence_from_refseq_entry(entry)

            # search for a match in gene list
            gene_found = False
            isoform_processed = False
            for gene in list_of_gene_objects:
                if isoform_processed:
                    break
                if HCGN_found:
                    if gene.HGNC == HGNC_ID:
                        gene_found = True
                if NCBI_ID_found:
                    if gene.refseq_gene_ID == NCBI_ID:
                        gene_found = True
                if gene_found:
                    if isoform_processed:
                        break
                    # print(gene.HGNC, HGNC_ID)
                    if NP:
                        if type(gene.protein_sequence_isoform_collection) == list:
                            for isoform in gene.protein_sequence_isoform_collection:
                                if isoform_processed:
                                    break
                                # print(isoform.refseq_NP, NP_ID)
                                if isoform.refseq_NP == NP_ID:
                                    if isoform.protein_sequence == protein_sequence:
                                        isoform.refseq_NP_version = NP_version
                                        isoform.refseq_NM_version = NM_ID_version
                                        print('ID versions added')
                                        NP_ID_version_added += 1
                                        isoform_processed = True
                                        break
                                    else:
                                        print('same NP ID but not same sequence')
                                        NP_ID_not_same_sequence += 1
                                elif isoform.protein_sequence == protein_sequence:
                                    isoform.refseq_NP = NP_ID
                                    isoform.refseq_NP_version = NP_version
                                    isoform.refseq_NM_version = NM_ID_version
                                    same_sequence_added_IDs += 1
                                    isoform_processed = True

                            if not isoform_processed:
                                # print(gene.HGNC, HGNC_ID, NCBI_ID,gene.refseq_gene_ID,NP_ID)
                                gene.protein_sequence_isoform_collection.append(
                                    Protein_isoform(protein_sequence, gene_name=gene.ensembl_gene_symbol,
                                                    refseq_NM_version=NM_ID_version, refseq_NP=NP_ID,
                                                    refseq_NP_version=NP_version))
                                isoform_processed = True
                                # print('added isoform')
                                sequences_added += 1
                    elif XP:
                        isoform_processed = True
                        pass
                        # XP_added += 1
                        # isoform_processed = True
                        # if type(gene.protein_sequence_isoform_collection) == list:
                        #    gene.protein_sequence_isoform_collection.append(
                        #        Protein_isoform(protein_sequence,gene_name=gene.ensembl_gene_symbol,
                        #                        refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,
                        #                        refseq_XP_version=XP_version))
                        # else:
                        #    gene.protein_sequence_isoform_collection = [Protein_isoform(protein_sequence,gene_name=gene.ensembl_gene_symbol,
                        #                        refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,
                        #                        refseq_XP_version=XP_version)]
                        #    print('new collection XP')
                    elif YP:
                        isoform_processed = True
                        #YP_added += 1
                        #isoform_processed = True
                        #if type(gene.protein_sequence_isoform_collection) == list:
                        #    gene.protein_sequence_isoform_collection.append(
                        #        Protein_isoform(protein_sequence, gene_name=gene.ensembl_gene_symbol,
                        #                        refseq_YP_version=YP_version, refseq_YP=YP_ID,
                        #                        refseq_NC_version=NC_ID_version))
                        #else:
                        #    gene.protein_sequence_isoform_collection = [
                        #        Protein_isoform(protein_sequence, gene_name=gene.ensembl_gene_symbol,
                        #                        refseq_YP_version=YP_version, refseq_YP=YP_ID,
                        #                        refseq_NC_version=NC_ID_version)]
                        #    print('new collection YP')

            # no match in list of gene objects
            # make a new gene object
            if not gene_found:
                isoform_processed = True
                no_match += 1
                if NP:
                    if HCGN_found:
                        list_of_gene_objects.append(Gene(HGNC=HGNC_ID, refseq_gene_ID=NCBI_ID,
                                                         protein_sequence_isoform_collection=[
                                                             Protein_isoform(protein_sequence,
                                                                             refseq_NM_version=NM_ID_version,
                                                                             refseq_NP=NP_ID,
                                                                             refseq_NP_version=NP_version)]))
                    else:
                        list_of_gene_objects.append(Gene(refseq_gene_ID=NCBI_ID,
                                                         protein_sequence_isoform_collection=[
                                                             Protein_isoform(protein_sequence,
                                                                             refseq_NM_version=NM_ID_version,
                                                                             refseq_NP=NP_ID,
                                                                             refseq_NP_version=NP_version)]))

                elif XP:
                    pass
                    isoform_processed = True
                    # if HCGN_found:
                    #    list_of_gene_objects.append(Gene(HGNC=HGNC_ID,refseq_gene_ID=NCBI_ID,
                    #                                     protein_sequence_isoform_collection=[Protein_isoform(protein_sequence,refseq_XM_version=XM_ID_version,refseq_XP=XP_ID,refseq_XP_version=XP_version)]))
                    # else:
                    #    list_of_gene_objects.append(Gene(refseq_gene_ID=NCBI_ID,
                    #                                     protein_sequence_isoform_collection=[
                    #                                         Protein_isoform(protein_sequence, refseq_XM_version=XM_ID_version, refseq_XP=XP_ID,refseq_XP_version=XP_version)]))
                elif YP:
                    pass
                    isoform_processed = True
                    #if HCGN_found:
                    #    list_of_gene_objects.append(Gene(HGNC=HGNC_ID, refseq_gene_ID=NCBI_ID,
                    #                                     protein_sequence_isoform_collection=[
                    #                                         Protein_isoform(protein_sequence, refseq_YP=YP_ID,
                    #                                                         refseq_NC_version=NC_ID_version,
                    #                                                         refseq_YP_version=YP_version)]))
                    #else:
                    #    list_of_gene_objects.append(Gene(refseq_gene_ID=NCBI_ID,
                    #                                     protein_sequence_isoform_collection=[
                    #                                         Protein_isoform(protein_sequence, refseq_YP=YP_ID,
                    #                                                         refseq_NC_version=NC_ID_version,
                    #                                                         refseq_YP_version=YP_version)]))

        print('total entries: ', len(splittext))
        print('not any type: ', not_any_type)
        print('no matches: ', no_match)
        print('NP sequences added: ', sequences_added)
        print('match but no isoforms:', match_but_no_isoforms)
        print('HCGN count', HCGN_count)
        print('NCBI count', NCBI_count)
        print('NP_ID_not_same_Seq', NP_ID_not_same_sequence)
        print('NP_Version added', NP_ID_version_added)
        print('same sequence found, added IDs', same_sequence_added_IDs)
        print('XP added', XP_added)
        print('YP added', YP_added)


    @staticmethod
    def extract_protein_sequence_from_refseq_entry(entry):
        '''extract protein sequence and format it correctly'''
        section = re.split("ORIGIN", entry)[1]
        without_newline = re.sub('\n', '', section)
        without_newline_and_whitespice = re.sub('\s', '', without_newline)
        without_numbers = re.sub("\d+", '', without_newline_and_whitespice)
        sequence_of_AA_acronym = '[A-Z]+'
        minimal_length_of_AA_seq = '[A-Z]{10}'
        raw_AA_seq_list = re.findall(sequence_of_AA_acronym + minimal_length_of_AA_seq, without_numbers.upper())
        if len(raw_AA_seq_list) >= 1:
            return raw_AA_seq_list[0]  # string
        else:
            return 'no AA sequence found'

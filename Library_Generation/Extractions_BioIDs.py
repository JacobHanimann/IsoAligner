import re


class Get_Bio_ID():
    pass


    @staticmethod
    def get_bio_IDs_with_regex(ID_type, string):
        'generic functions to extract certain ID types from different databases'
        version = False
        # Ensembl
        if ID_type == 'ensembl_ensg':
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
        elif ID_type == 'ensembl_ense':
            pattern = 'ENSE\d{11}'
        elif ID_type == 'ensembl_ense_version':
            pattern = 'ENSE\d+\.\d+'
            version = True

        # Refseq
        # known
        elif ID_type == 'refseq_NM':
            pattern = 'NM_\d+'
        elif ID_type == 'refseq_rna_version':
            pattern = 'NM_\d+\.\d+'
            version = True
        elif ID_type == 'refseq_prot':
            pattern = 'NP_\d+'
        elif ID_type == 'refseq_prot_version':
            version = True
            pattern = 'NP_\d+\.\d+'
        elif ID_type == 'refseq_prot_predict':
            # model refseq
            pattern = 'XP_\d+'
        elif ID_type == 'refseq_prot_predict_version':
            pattern = 'XP_\d+\.\d+'
            version = True
        elif ID_type == 'refseq_rna_predict':
            pattern = 'XM_\d+'
        elif ID_type == 'refseq_rna_predict_version':
            pattern = 'XM_\d+\.\d+'
            version = True
        # mitchondrium
        elif ID_type == 'refseq_mitocho':
            pattern = 'YP_\d+'
        elif ID_type == 'refseq_mitocho_version':
            version = True
            pattern = 'YP_\d+\.\d+'
        # refseq chromosome
        elif ID_type == 'refseq_chromosome':
            pattern = 'NC_\d+\.\d+'
        elif ID_type == 'refseq_chromosome_version':
            version = True
            pattern = 'NC_\d+\.\d+'

        # Uniprot IDs
        elif ID_type == 'uniprot_accession':
            pattern = '\|[OPQ][0-9][0-9A-Z]{3}[0-9]\||\|[A-NR-Z][0-9][A-Z][A-Z,0-9]{2}[0-9]\||\|[A-N,R-Z][0-9][A-Z][A-Z,0-9]{2}[0-9][A-Z][A-Z,0-9]{2}[0-9]\|'
        elif ID_type == 'uniprot_uniparc':
            pattern = 'UPI[0-9A-F]+'

        if ID_type == "gene_name":
            pattern = "\|[^\|\n]+\n"
            match_list = re.findall(pattern, string)
            if not match_list:  # if list is empty
                return 'not found'
            else:
                # remove \n with regex
                return match_list[0][1:-1]  # remove \n

        # HGNC
        elif ID_type == 'HGNC':
            pattern = "HGNC:\d+"

        # execute regular expression
        match_list = re.findall(pattern, string)
        if not match_list:  # if list is empty
            return 'not found'
        elif len(match_list) == 1:
            if ID_type == "uniprot_accession":
                return match_list[0][1:-1]
            else:
                return match_list[0]
        else:
            if version == False:
                return match_list[0]
            else:
                return match_list

class Exon():
    def __init__(self,exon_start_end=None,exon_length=None,ENSE=None, ENST =None, exon_number=None):
        self.exon_start_end = exon_start_end
        self.exon_length_in_AA = exon_length
        self.ENSE = ENSE
        self.ENST = ENST
        self.exon_number = exon_number


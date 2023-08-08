from Bio import Align


MODE = 'local'
MATCH_SCORE = 1.0
MISMATCH_SCORE = -2.0
GAP_SCORE = -5


class LocalAlignment(Align.PairwiseAligner):
    """PairwiseAligner pour alignement local."""
    def __init__(self, mode=MODE, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE, gap_score=GAP_SCORE):
        super().__init__() # Calls the __init__ method of the parent class (PairwiseAligner) good practice to do it first before the child class __init__ is performed
        self.mode = mode # Define default attributes or data
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

    def run(self, template, query):
        """
        Démarre un alignement local entre deux objets SeqRecord à l'aide de l'objet LocalAlignment.

        Args:
            template: (objet: Seq) qui sert de template pour l'alignement.
            query: (objet: Seq) est la séquence à aligner.

        Returns:
            Alignement de deux séquences sous la forme d'un objet Bio.Align.PairwiseAlignment.

        """
        score = self.score(template, query)
        score_reverse = self.score(template, query.reverse_complement())
        if score > score_reverse:
            alignment = super().align(template, query)[0]
        else:
            alignment = super().align(template, query.reverse_complement())[0]
        return alignment

    def get_strand(self, template, query):
        """Renvoie le sens du read."""
        score = self.score(template, query)
        score_reverse = self.score(template, query.reverse_complement())
        if score > score_reverse:
            return "+"
        else:
            return "-"

"""
distance.py

holds distance classes

"""

import editdistance  # levenshtein distance
import numpy


# abstract base class for distances
class DistMat(object):

    """
    distances
    """

    def __init__(self, sequence=None, length=None, fastaid=None, folds=None, subopt_folds=None):
        self.sequence = sequence
        self.length = length
        self.fastaid = fastaid
        self.subopt_folds = subopt_folds


class LevDistance(DistMat):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # for instantiating from DotBracket class
    def from_DotBracket(self, DBclass):
        folded_strucs = DBclass.subopt_folds
        # remove counter and sequence keys from subopt dict
        folded_strucs.pop('counter', None)
        folded_strucs.pop('sequence', None)

        distance_mat = np.zeros((DBclass.length, DBclass.length))

    @staticmethod
    def compute_Ldist():
        pass

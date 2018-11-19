"""
distance.py

holds distance classes

"""

import editdistance  # levenshtein distance
import numpy as np
import time

# abstract base class for distances


class DistMat(object):

    """
    distances
    """

    def __init__(self, sequence=None,
                 length=None,
                 fastaid=None,
                 folds=None,
                 subopt_folds=None,
                 distances=None):
        self.sequence = sequence
        self.length = length
        self.fastaid = fastaid
        self.subopt_folds = subopt_folds
        self.distances = distances


class LevDistance(DistMat):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # for instantiating from DotBracket class

    @classmethod
    def from_DotBracket(cls, DBclass):
        folded_strucs = DBclass.subopt_folds
        folded_strucs = list(folded_strucs.values())

        # number of folded structures
        num_folded_structs = len(folded_strucs)

        # initialize distance array
        distance_mat = np.zeros((num_folded_structs, num_folded_structs))

        print("computing distances..")
        ts = time.time()
        for indx_i in range(len(folded_strucs)):
            for indx_j in range(indx_i, len(folded_strucs)):

                # compute distance
                ldist = editdistance.eval(folded_strucs[indx_i], folded_strucs[indx_j])

                # save it to array
                distance_mat[indx_i][indx_j] = ldist
        distance_mat += distance_mat.T

        te = time.time()
        print("finished computing {} levenshtein distances in {:.3f} sec".format(np.count_nonzero(distance_mat) + num_folded_structs,
                                                                                 te - ts))

        return cls(sequence=DBclass.sequence, length=DBclass.length,
                   fastaid=DBclass.fastaid, distances=distance_mat, subopt_folds=DBclass.subopt_folds)

    @staticmethod
    def compute_Ldist(folded_strucs):
        """takes in a dictionary of N folded structures

        Args:
            folded_strucs (dict): dictionary of N folded structures

        Returns:
            distance_mat (numpy array): NxN numpy array of distances
        """
        folded_strucs = list(folded_strucs.values())
        # number of folded structures
        num_folded_structs = len(folded_strucs)

        # initialize distance array
        distance_mat = np.zeros((num_folded_structs, num_folded_structs))

        print("computing distances..")
        ts = time.time()
        for indx_i in range(len(folded_strucs)):
            for indx_j in range(indx_i, len(folded_strucs)):

                # compute distance
                ldist = editdistance.eval(folded_strucs[indx_i], folded_strucs[indx_j])

                # save it to array
                distance_mat[indx_i][indx_j] = ldist

        te = time.time()
        print("finished computing {} levenshtein distances in {:.3f}".format(np.count_nonzero(distance_mat) + num_folded_structs,
                                                                             te - ts))

        return distance_mat

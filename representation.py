"""
Representations.py

holds classes for different ways of representing the secondary structures
"""


class Representation(object):

    """
    representation
    """

    def __init__(self, sequence=None, length=None, fastaid=None, folds=None, subopt_folds=None):
        self.sequence = sequence
        self.length = length
        self.fastaid = fastaid
        self.subopt_folds = subopt_folds


# representation child class that will convert self data into
# other types of data

class DotBracket(Representation):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_structure(cls, folded_struc):
        # unpack feed information into attributes
        cls.sequence = folded_struc.sequence
        cls.fastaid = folded_struc.fastaid
        if folded_struc.subopt_folds is None:
            raise ValueError("Error: Not folds found in {}".format(folded_struc))
        else:
            cls.subopt_folds = folded_struc.subopt_folds
            cls.length = folded_struc.subopt_folds['counter']
        return cls

    def convert_to_seq():
        pass

"""
Secondary Structure classes 

Going to call file: "structure.py"
RNA folding

To do:
- Need to automatically update length attribute when a sequence attribute is set
"""


import RNA  # RNAfold python wrapper
#import abc

# abstract class for any 2D structure


class SecondaryStructure(object):

    """
    class for the secondary structure

    Attributes:
                    fastid (int): Description
                    length (int): length of primary monomer sequence
                    sequence (str): monomer sequence
    """

    def __init__(self, sequence=None, length=None, fastaid=None):
        self.sequence = sequence
        self.fastaid = fastaid

        if sequence is not None:
            self.length = len(sequence)
        else:
            self.length = length


class RNAfolded(SecondaryStructure):

    def __init__(self, subopt_folds=None, *args, **kwargs):
        self.subopt_folds = subopt_folds
        super().__init__(*args, **kwargs)

    def __str__(self):
        if self.subopt_folds == None:
            num_folds = 0
        else:
            num_folds = len(self.subopt_folds)

        return "{}, sequence:{}, length:{}, folds:{}".format(self.__class__,
                                                             self.sequence,
                                                             self.length,
                                                             num_folds)

    # helper function for RNAfold
    @classmethod
    def folding_to_dict(cls, structure, energy, data):
        """
        helper function for the fold_from_str() function
        """
        # dictionary to hold the folding
        # subopt_folds = {'counter': 1, 'sequence': sequence}

        if not structure == None:
            # struct_id
            stru_id = "structure_{}_{:.4f}".format(data['counter'], energy)

            # save structure
            data[stru_id] = structure

            # increase structure counter
            data['counter'] = data['counter'] + 1

    def fold(self, ewindow=500):
        if self.sequence == None:
            print("Error: No sequence found!")
        else:
            # Set global switch for unique ML decomposition
            RNA.cvar.uniq_ML = 1

            # initalize dictionary to store folds
            self.subopt_folds = {'counter': 1, 'sequence': self.sequence}

            # Create a 'fold_compound' for our sequence
            RNAfold_compound = RNA.fold_compound(self.sequence)

            print('folding RNA...')
            # energy window is in dacal units (500 dacal/mol = 5 kcal/mol)
            RNAfold_compound.subopt_cb(ewindow,
                                       self.folding_to_dict,
                                       self.subopt_folds)
            print('finished folding!')

            # returns dictionary of folds
            return self.subopt_folds

    # for instatiating the class from a Sequence class
    @classmethod
    def fold_from_str(cls, str_seq, ewindow=500):
        """Folding from a sequence string """

        subopt_folds = {'counter': 1, 'sequence': str_seq}

        # Set global switch for unique ML decomposition
        RNA.cvar.uniq_ML = 1

        # create fold compound for our sequence
        RNAfold_compound = RNA.fold_compound(str_seq)

        # fold the sequence
        # energy window is in dacal units (500 dacal/mol = 5 kcal/mol)

        RNAfold_compound.subopt_cb(ewindow,
                                   cls.folding_to_dict,
                                   subopt_folds)

        return cls(sequence=str_seq, length=len(str_seq), subopt_folds=subopt_folds)

"""
Classes for the primary structure

Going to call file: "sequence.py"

To do:
- extend the RNA child class for RNA-specific functionality
- make PrimeStruc an abstract class
"""


# abstract class for any 1D sequence

class PrimStruc(object):

    """
    class for the intial data

    Attributes:
        composition (dict): statistics of sequence
        length (int): count of monomers in sequence
        sequence (str): sequence as string
    """

    def __init__(self,
                 sequence=None,
                 length=None,
                 fastaid=None):
        self.sequence = sequence
        self.length = length
        self.composition = {}
        # self.fastaid = fastaid

    def comp_stats(self, percent=False):
        """
        statistics from the sequence
        """

        # check that sequence is there
        if len(self.sequence) == 0:
            raise ValueError("ERROR: no sequence found")

        else:
            print("getting composition stats")
            # iterate over characters and save them
            for char in self.sequence:
                self.composition[char] = self.sequence.count(char)

            # check percent bool
            if percent:
                for char, count in self.composition.items():
                    self.composition[char] = count / self.length

        return self.composition

    # @classmethod
    def from_txt(self, filename):
        """import sequence from txt file
        s
        Args:
            filename (txt): text file with only sequence inside

        Returns:
            sequence (str): attribute of class
        """
        with open(filename) as f:
            for line in f:
                seq = line
                if 'str' in line:
                    break
        self.sequence = seq
        self.length = len(seq)

        return self.sequence

    @classmethod
    def from_str(self, seqstring):
        """import class from str

        gives a new instantiaion whith sequence provided
        Args:
            seqstring (str): sequence as seqstring
        """
        primarystructure_inst = PrimStruc(sequence=seqstring,
                                          length=len(seqstring))

        return primarystructure_inst


class RNAseq(PrimStruc):

    """
    RNA sequence 
    """

    # def __int__(self, sequence, length, fastaid):
    #    super().__init__(*args, *kwargs)

    pass

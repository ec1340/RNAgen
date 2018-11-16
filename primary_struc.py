"""
Classes for the primary structure
"""

# abstract class for any 1D sequence
class PrimStruc():

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

    @classmethod
    def from_txt(self, filename):
        """import sequence from txt file
        
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
    def from_str(self, string):
        """import sequence from str

        Args:
            string (str): sequence as string
        """
        primarystructure_inst = PrimStruc(sequence=string,
        									length = len(string))

        return primarystructure_inst


# class RNA_seq()

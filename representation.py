"""
Representations.py

holds classes for different ways of representing the secondary structures
"""
import numpy as np
import re
import scipy


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

    def __init__(self, energies=None, *args, **kwargs):
        self.energies = energies
        super().__init__(*args, **kwargs)

    @classmethod
    def from_structure(cls, folded_struc):
        # unpack feed information into attributes
        cls.sequence = folded_struc.sequence
        cls.fastaid = folded_struc.fastaid
        cls.length = folded_struc.length

        # subopt_folds attribute holds dot bracket information

        if folded_struc.subopt_folds is None:
            raise ValueError("Error: No folds found in {}".format(folded_struc))
        else:
            cls.subopt_folds = folded_struc.subopt_folds
        return cls

    @staticmethod
    def expander(i_start, j_end, length):
        """Expands short notation of contacts

        helper function to create (i,j) tuples where bases i and j make contact
        in the RNA secondary structure  - used in get_contacts_from_txt()


        output: list of contact (i,j) tuples

        Args:
            i_start (int): starting index
            j_end (int): ending index
            length (int): length of pairing stem

        Returns:
            pairs (dict): dictionary of pairings
        """

        # create list of paired i indices
        i_s = list(np.arange(i_start, i_start + length))

        # create list of corresponding j indices:
        j_s = list(np.arange(j_end, j_end + length))

        # form list of tuple pairs
        pairs = list(zip(i_s, j_s))

        return pairs

    @classmethod
    def get_contacts_from_txt(cls, filename):
        """Reads txt file of RNA contacts into dictionary

        reads RNA secondary structure contact file into a dictionary
        where each key is a structure tag and each value is a list of tuples

        Args:
            filename (txt file): file

        Returns:
            struc_contact_dict (dict): dictionary of pairings
            struct_Es (list):  M length list of energies

        """
        struc_contact_dict = {}
        contact_list = None
        current_key = 0

        struct_Es = []

        # regex pattern
        pattern = "(Structure)\s+(\d+)\s+(-\d+\.\d+)\s+(\d+\.\d+)"
        regex = re.compile(pattern)

        with open(filename) as f:
            for ind, line in enumerate(f):
                m = re.search(string=line, pattern=regex)

                # structure heading lines
                if m is not None:

                    # save list to current key before we re-initialize the list
                    # for the next key
                    if contact_list is not None:
                        struc_contact_dict[current_key] = contact_list

                    contact_list = []

                    # create new key
                    current_key = "_".join(str(x) for x in list(m.groups()))

                    # add energy to list
                    struct_Es.append(float(m.groups()[2]))

                # contact indices lines
                if m is None:

                    # transform line to list of integers that will go into expander fcn
                    expander_input = [int(i) for i in line.split()]

                    contacts = cls.expander(expander_input[0], expander_input[1], expander_input[2])

                    # add contacts to list for this structure key
                    contact_list = contact_list + contacts

                # checks for end of file
                if 'str' in line:
                    break

            # add last entry
            struc_contact_dict[current_key] = contact_list

        return struc_contact_dict, struct_Es

    @classmethod
    def from_txt(cls, filename, seq_length):
        """reads file into dotbracket notation

        Args:
            filename (TYPE): Description
        """
        contact_dict, energies = cls.get_contacts_from_txt(filename)

        # dot bracket list for all structurs
        all_struct_dot_bracket_dict = {}

        # iteratre over len and keys of the contact dict
        for struc_indx, struc_key in enumerate(contact_dict.keys()):

            # current structure in the contact dict
            current_entry = contact_dict[struc_key]

            # dot bracket list for current entry - # dot initialized
            current_entry_DB_list = "." * seq_length
            current_entry_DB_list = list(current_entry_DB_list)

            # iterate over the pair tuples in the
            for pair in current_entry:

                # get indices from tuple
                open_par_indx = int(pair[0])
                clos_par_indx = int(pair[1])

                # open par
                current_entry_DB_list[open_par_indx] = "("

                # closed par
                current_entry_DB_list[clos_par_indx] = "("

            # save to dot bracket dict
            all_struct_dot_bracket_dict[struc_key] = ''.join(current_entry_DB_list)

        return cls(subopt_folds=all_struct_dot_bracket_dict, length=seq_length)


class OneChannel(Representation):

    def __init__(self, energies=None, *args, **kwargs):
        self.energies = energies
        super().__init__(*args, **kwargs)

    @staticmethod
    def flatten_mat(contact_mat):
        """flattens the contact matrix by taking the

        Args:
            contact_mat (MxNxN array): NxN contact matrices
            binary (bool, optional): converts all values >1 to 1

        Returns:
            flat_contact_mat: MxNx1 array
        """

        colsum_mat = 0
        for mat in contact_mat:
            # check for first matrix
            if type(colsum_mat) != np.ndarray:
                # flatten the matrix down by getting column sums
                colsum_mat = np.sum(mat, axis=0)

            # for any matrix except the first one
            else:
                colsum_mat = np.concatenate((colsum_mat,
                                             np.sum(mat, axis=0)),
                                            axis=0)

        flat_contact_mat = colsum_mat.reshape(contact_mat.shape[0], -1, 1)

        for indx, entry in enumerate(flat_contact_mat):
            # convert non-zero entries to 1
            flat_contact_mat[indx] = np.where(entry != 0, 1, 0)

        return flat_contact_mat

    @staticmethod
    def contact_dict_to_matrix(contact_dict):
        """Creates indicator matrix from contacts dict

        takes in a dictionary of tuples (i, j) contacts and
        writes them into corresponding indicator values in
        an NxN matrix, where N=len of RNA sequence


        Args:
            contact_dict (dictionary): dictionary of pairings

        Returns:
            mat_array: returns an MxN matrix
        """

        # initialize np array
        mat_array = []

        for struc_id, stru_contacts in contact_dict.items():
            i = 0

            # initalize current contact matrix
            current_array = np.zeros([500, 500])

            # read contact tuples into array
            for cont in stru_contacts:
                i, j = cont
                current_array[i][j] = 1

            # add this matrix to array of matrices
            mat_array.append(current_array)

            # move to next matrix (structure)
            i += 1

        mat_array = np.array(mat_array)

        return mat_array

    @classmethod
    def from_txt(cls, filename, seq_length):
        """Reads txt file of RNA contacts into dictionary

        reads RNA secondary structure contact file into a dictionary
        where each key is a structure tag and each value is a list of tuples

        Args:
            filename (txt file): file

        Returns:
            struc_contact_dict (dict): dictionary of pairings
            struct_Es (list):  M length list of energies

        """
        struc_contact_dict = {}
        contact_list = None
        current_key = 0

        struct_Es = []

        # regex pattern
        pattern = "(Structure)\s+(\d+)\s+(-\d+\.\d+)\s+(\d+\.\d+)"
        regex = re.compile(pattern)

        with open(filename) as f:
            for ind, line in enumerate(f):
                m = re.search(string=line, pattern=regex)

                # structure heading lines
                if m is not None:

                    # save list to current key before we re-initialize the list
                    # for the next key
                    if contact_list is not None:
                        struc_contact_dict[current_key] = contact_list

                    contact_list = []

                    # create new key
                    current_key = "_".join(str(x) for x in list(m.groups()))

                    # add energy to list
                    struct_Es.append(float(m.groups()[2]))

                # contact indices lines
                if m is None:

                    # transform line to list of integers that will go into expander fcn
                    expander_input = [int(i) for i in line.split()]

                    contacts = DotBracket.expander(expander_input[0], expander_input[1], expander_input[2])

                    # add contacts to list for this structure key
                    contact_list = contact_list + contacts

                # checks for end of file
                if 'str' in line:
                    break

            # add last entry
            struc_contact_dict[current_key] = contact_list

        # convert dictionary to matrix
        contact_mat = cls.contact_dict_to_matrix(struc_contact_dict)

        # convert matrix to flat matrix
        flat_contact_mat = cls.flatten_mat(contact_mat)

        return cls(subopt_folds=flat_contact_mat, length=seq_length)

    def smoothen(self):

        return SmoothOneChannel(subopt_folds="pie")


class SmoothOneChannel(OneChannel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

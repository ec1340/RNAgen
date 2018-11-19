"""
Examples of using the RNAfold python interfaces

examples are taken from:
https://github.com/ViennaRNA/ViennaRNA/tree/master/examples/Python

"""


import RNA

# # -------------------
# # Example 1a
# # -------------------

# # The RNA sequence
# seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"

# # compute minimum free energy (MFE) and corresponding structure
# (ss, mfe) = RNA.fold(seq)

# # print output
# print("%s\n%s [ %6.2f ]" % (seq, ss, mfe))


# # -------------------
# # Example 2
# # -------------------

# # The RNA sequence alignment
# sequences = [
#     "CUGCCUCACAACGUUUGUGCCUCAGUUACCCGUAGAUGUAGUGAGGGU",
#     "CUGCCUCACAACAUUUGUGCCUCAGUUACUCAUAGAUGUAGUGAGGGU",
#     "---CUCGACACCACU---GCCUCGGUUACCCAUCGGUGCAGUGCGGGU"
# ]

# # compute the consensus sequence
# cons = RNA.consensus(sequences)

# # predict Minmum Free Energy and corresponding secondary structure
# (ss, mfe) = RNA.alifold(sequences)

# # print output
# print("%s\n%s [ %6.2f ]" % (cons, ss, mfe))

# # -------------------
# # Example 3
# # -------------------


# # The RNA sequence
# seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"

# # create a new model details structure
# md = RNA.md()

# # change temperature and dangle model
# md.temperature = 20.0  # 20 Deg Celcius
# md.dangles = 1    # Dangle Model 1

# # create a fold compound
# fc = RNA.fold_compound(seq, md)

# # predict Minmum Free Energy and corresponding secondary structure
# (ss, mfe) = fc.mfe()

# # print sequence, structure and MFE
# print("%s\n%s [ %6.2f ]\n" % (seq, ss, mfe))

# # -------------------
# # Example 4
# # -------------------

# """

# seq1 = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCAA"

# # Turn-off dangles globally
# RNA.cvar.dangles = 0

# # Data structure that will be passed to our MaximumMatching() callback with two components:
# # 1. a 'dummy' fold_compound to evaluate loop energies w/o constraints, 2. a fresh set of energy parameters
# mm_data = { 'dummy': RNA.fold_compound(seq1), 'params': RNA.param() }

# # Nearest Neighbor Parameter reversal functions
# revert_NN = {
#     RNA.DECOMP_PAIR_HP:       lambda i, j, k, l, f, p: - f.eval_hp_loop(i, j) - 100,
#     RNA.DECOMP_PAIR_IL:       lambda i, j, k, l, f, p: - f.eval_int_loop(i, j, k, l) - 100,
#     RNA.DECOMP_PAIR_ML:       lambda i, j, k, l, f, p: - p.MLclosing - p.MLintern[0] - (j - i - k + l - 2) * p.MLbase - 100,
#     RNA.DECOMP_ML_ML_STEM:    lambda i, j, k, l, f, p: - p.MLintern[0] - (l - k - 1) * p.MLbase,
#     RNA.DECOMP_ML_STEM:       lambda i, j, k, l, f, p: - p.MLintern[0] - (j - i - k + l) * p.MLbase,
#     RNA.DECOMP_ML_ML:         lambda i, j, k, l, f, p: - (j - i - k + l) * p.MLbase,
#     RNA.DECOMP_ML_UP:         lambda i, j, k, l, f, p: - (j - i + 1) * p.MLbase,
#     RNA.DECOMP_EXT_STEM:      lambda i, j, k, l, f, p: - f.E_ext_loop(k, l),
#     RNA.DECOMP_EXT_STEM_EXT:  lambda i, j, k, l, f, p: - f.E_ext_loop(i, k),
#     RNA.DECOMP_EXT_EXT_STEM:  lambda i, j, k, l, f, p: - f.E_ext_loop(l, j),
#     RNA.DECOMP_EXT_EXT_STEM1: lambda i, j, k, l, f, p: - f.E_ext_loop(l, j-1),
#             }

# # Maximum Matching callback function (will be called by RNAlib in each decomposition step)
# def MaximumMatching(i, j, k, l, d, data):
#     return revert_NN[d](i, j, k, l, data['dummy'], data['params'])

# # Create a 'fold_compound' for our sequence
# fc = RNA.fold_compound(seq1)

# # Add maximum matching soft-constraints
# fc.sc_add_f(MaximumMatching)
# fc.sc_add_data(mm_data, None)

# # Call MFE algorithm
# (s, mm) = fc.mfe()

# # print result
# print("%s\n%s (MM: %d)\n" %  (seq1, s, -mm))
# """
# print("not running since it kills the jupyter kernel")

# # -------------------
# # Example 5
# # -------------------


# sequence = "CGCAGGGAUACCCGCG"

# # create new fold_compound object
# fc = RNA.fold_compound(sequence)

# # compute minimum free energy (mfe) and corresponding structure
# (ss, mfe) = fc.mfe()

# # print output
# print("%s [ %6.2f ]" % (ss, mfe))


# -------------------
# Example 6
# -------------------

sequence = "GGGGAAAACCGACGACGACAGCAGGCAGGCAGCGAGGCAGGCAGGCGCCCGCGGAACCC"

# Set global switch for unique ML decomposition
RNA.cvar.uniq_ML = 1

subopt_data = {'counter': 1, 'sequence': sequence}

# Print a subopt result as FASTA record


def print_subopt_result(structure, energy, data):
    if not structure == None:
        # print(">subopt %d" % data['counter'])
        # print("%s" % data['sequence'])
        # print("%s [%6.2f]" % (structure, energy))

        # struct_id
        stru_id = "structure-{}-{:.4-f}".format(data['counter'], energy)
        # stru_id = "structure-" + str(data['counter']) + str(energy)

        # save structure
        subopt_data[stru_id] = structure

        # increase structure counter
        data['counter'] = data['counter'] + 1



# Create a 'fold_compound' for our sequence
a = RNA.fold_compound(sequence)

# Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
# the MFE and print each structure using the function above
b = a.subopt_cb(100, print_subopt_result, subopt_data)

print(subopt_data)

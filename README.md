# RNAgen


- using RNA_env conda environment

### Workflow

1) Parse sequence from file

2) Generate secondary structure (RNAfold)

3) Transform secondary structures to dot bracket notation

- should also include options for other representation methods

4) Use L distances between the secondary structures generated

5) Using the produced distance matrix, run PHATE or diffusion maps 

### Requirements

- Biopython

- RNAfold
-- installation instructions at  https://github.com/ViennaRNA/ViennaRNA
-- Route taken: 

#### RNA fold installation steps
- check that GNU compliler tools are available https://superuser.com/questions/383580/how-to-install-autoconf-automake-and-related-tools-on-mac-os-x-from-source
- check that activate developer path https://apple.stackexchange.com/questions/254380/macos-mojave-invalid-active-developer-path
- ./configure option for python 3: 
		./configure --disable-openmp --without-perl --without-python
- run 
		> make 
- run 
		> make install
- if RNA isn't found when you run >import RNA, then go to your .bash-profile and then add the location of the python 3 site packages to your PYTHONPATH variable https://leemendelowitz.github.io/blog/how-does-python-find-packages.html


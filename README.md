# RNAgen

### Quick Start (in progress)


### Requirements

- RNAfold
- numpy

#### RNA fold installation steps
- check that GNU compliler tools are installed:
	- https://superuser.com/questions/383580/how-to-install-autoconf-automake-and-related-tools-on-mac-os-x-from-source
- check that xcode developer tools are installed
	- if not an "inactive developer path error" will be raised (fix if recently ugraded to OS X Mojave: https://apple.stackexchange.com/questions/254380/macos-mojave-invalid-active-developer-path)

- follow quick start installation instructions at https://github.com/ViennaRNA/ViennaRNA 
- NOTE: ./configure option for python 3: 
		./configure --disable-openmp --without-perl --without-python
- if RNA isn't found when you run >import RNA, then go to your .bash-profile and then add the location of the python 3 site packages to your PYTHONPATH variable https://leemendelowitz.github.io/blog/how-does-python-find-packages.html


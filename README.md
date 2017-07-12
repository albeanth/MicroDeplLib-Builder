## MAMMOTH Depletion Library Builder
This script will build depletion libraries for the MAMMOTH reactor physics code package. It reads in the following [ENDF7.1] (http://www.nndc.bnl.gov/endf/b7.1/download.html) sublibraries:

- Neutron Reaction Sublibrary
- Decay Reaction Sublibrary
- Spontaneous Fission Product Yields Reaction Sublibrary
- Neutron Induced Fission Product Yields Sublibrary

Python3.x+ is required. Originally built with Python 3.5.1 using a full Anaconda build on OSX 10.11.5. The following python packages are required:

- sys,os
- os.path
- linecache
- re
- numpy
- xml.etree.ElementTree
- zipfile
- urlopen

The package [colorama] (https://pypi.python.org/pypi/colorama) is used for pretty outputting. If not used, print statements will be standard.

#### Execution
The MicroDeplLib-Build.py script will initially look for the required NNDC ENDF7.1 files. If they do not exist, they a directory called **./ENDF7.1** will be created. The required sublibraries will be downloaded as zip files, unzipped, and moved into **./ENDF7.1**. Once moved, the original zip files are deleted. The execution of the program is as follows:

    > python MicroDeplLib-Build.py

Once it is complete, the following will be outputted *Done building library.*. To parse the created xml file into a human readable format, execute **convert.sh**:

    > ./convert.sh

The main python script will take on the order of 5 minutes to run (machine dependent). The outputted xml file after parsing will be 2.1MB.

#### Known Issues
*1. Spontaneous Fission*
There are more isotopes in the decay sublibrary that identify as capable of undergoing spontaneous fission  (SF) as there are in the SF sublibrary. This leads to several isotopes having unknown SF yields.

*2. Neutron Induced Fission*
As in the SF case, there are more isotopes in the neutron reaction sublibrary that identify as fissionable than there are in the neutron induced fission product yield library. This leads to several isotopes having unknown neutron induced fission yields.

*3. Decay Data*
There are 9 cases in which isotopes decay to unknown products. These cases are hard coded in the **TranslateDecayMode** function in MicroDeplLib-Build.py.

If any new bugs are found, please report them to Tony Alberti at albertia@oregonstate.edu

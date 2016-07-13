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

The package [colorama] (https://pypi.python.org/pypi/colorama) is used for pretty outputting. If not used, print statements will need to be adjusted. 

#### Execution
You will need to have the main MicroDeplLib-Build.py and translators.py files in the same directory (or at least having the MicroDeplLib-Build.py knowing where to find translators.py). Furthermore, you will need to adjust the pathing for the decay sublibraries within MicroDeplLib-Build.py. Once you have all of the proper paths set up, you can execute the program as follows.

    > python MicroDeplLib-Build.py
    > ./convert.sh

The python script will take on the order of 5 minutes to run (machine dependent). The convert.sh file is used to parse out the xml file and make it human readable. The outputted xml file after parsing will be 2.1MB. 

#### Known Issues
*1. Spontaneous Fission* 
There are more isotopes in the decay sublibrary that identify as capable of undergoing spontaneous fission  (SF) as there are in the SF sublibrary. This leads to several isotopes having unknown SF yields.

*2. Neutron Induced Fission*
As in the SF case, there are more isotopes in the neutron reaction sublibrary that identify as fissionable than there are in the neutron induced fission product yield library. This leads to several isotopes having unknown neutron induced fission yields. 

*3. Decay Data*
There are 9 cases in which isotopes decay to unknown products. These cases are hard coded in the **TranslateDecayMode** function in MicroDeplLib-Build.py. 

If any new bugs are found, please report them to Tony Alberti at albertia@oregonstate.edu

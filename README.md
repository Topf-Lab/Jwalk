# Jwalk

#===============================================================================
#     This file is part of Jwalk.
#     
#     Jwalk - A tool to calculate the solvent accessible surface distance (SASD) 
#     between crosslinked residues.
#     
#     Copyright 2016 Jwalk Inventor and Birkbeck College University of London.
#                          The Jwalk Inventor is: Josh Bullock
# 
# 
#     Jwalk is available under Public Licence.
#     This software is made available under GPL V3
#
#     Please cite your use of Jwalk in published work:
#     
#     J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#     The importance of non-accessible crosslinks and solvent accessible surface distance
#     in modelling proteins with restraints from crosslinking mass spectrometry. 
#     Molecular and Cellular Proteomics (15) pp.2491–2500
#
#===============================================================================

INSTALLATION:
Open a terminal session and a move to the directory Jwalk is in.

type "setup.py install"

RUNNING JWALK:

Jwalk currently only runs on Unix or Linux based operating systems (sorry !)

type "jwalk" plus any of the below flags:

  -h, --help        show this help message and exit
  -lys              calculate lysine crosslinks (default)
  -xl_list XL_LIST  calculate crosslinks from input list (see Examples)
  -i I              specify input pdb: -i <inputfile.pdb>   (default runs on every .pdb in directory)
  -aa1 AA1          specify starting amino-acid via three letter code eg. HIS
  -aa2 AA2          specify starting amino-acid via three letter code eg. TYR
  -surface          use higher accuracy method to calculate solvent accessibility - slower

if no flags are added, Jwalk will calculate the SASDs between all lysine residues on
every .pdb file in the directory.

CROSSLINK LIST INPUT:

Crosslinks should be listed in a pipe-delimited .txt file as so:

res1|chain1|res2|chain2|

eg. 34|A|56|B|

If your pdb file has no chain, insert a lowercase x

eg. 34|x|56|x|

See the example file for a sample list.

OUTPUT:

Two output files are generated in ./Jwalk_results/ 

<pdb>_crosslink_list.txt

<pdb>_crosslinks.pdb

In order to visualize the crosslinks in Chimera select Actions>Atoms/Bonds>show
Better representation in chimera in progress ...

Please try the test in JWalk/Examples/

Happy crosslinking !

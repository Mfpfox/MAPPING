------------------------------
------------------------------
------------------------------
------------------------------
#repo where i dowloaded origianl scripts:
# http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/

# origianl project dir : CODE_DATA/db_EXAC/KRASproject/

# DONT DO THIS 1. loadBfactors.py
1. data2bfactor.py
2. color_b.py

1. loadBfactors.py # used with CASP8, multi chains
2. alpha2all.py

---MAKING CADD SCORES FILE 
# in pymol set dir
[1] cd: now in '/Users/mariapalafox/Desktop/MAPPING/PDB_mapped_score'

# making atom specific file
[2] on PDB website, download the 'PDB Format (gz)', unzip, change extension to .txt instead of .pdb; .txt file, copy rows beginning with 'ATOM' into new file names '2BH9_ATOM.txt'

# NOTE- only include ATOM with corresponding positions that you are mapping CADD to, for KRAS, i only went to pos 172, trimmed at end
# NOTE- deleted row with TER instead of ATOM



# make column formating nicer for loading into pandas
[3] column -t atomfile.txt > formattedatom_file.txt


[4] load in pandas df
g6pdb = 'g6pd_2bh9_ATOM_start28pos.txt'
c8pdb = 'casp8_3kjn_ATOM_allpos.txt'
df6 = pd.read_csv(g6pdb, delimiter=r"\s+",header=None)
df8 = pd.read_csv(c8pdb, delimiter=r"\s+",header=None)
# formated column
df8.to_csv("RESIDUE_3kjn_casp8.csv", index=False)
df6.to_csv("RESIDUE_2bh9_g6pd.csv", index=False)
# had to delete some rows in casp8 that i think were there from binding small molecule


[5] add cadd average column to RESIDUE..csv file and check that the amino acids and positions are all correctly aligned
# NOTE using all positions in pdb ATOM file, CASP8 is perfect by g6pd pos 27 is V in pdb and H in ukb...that the only not matching pos tho. 


[6] create files to hardcode into data2bfactor.py
byResidue_2bh9_g6pd.txt
byResidue_3kjn_casp8.txt
------------------------------
------------------------------
------------------------------
------------------------------
------------------------------
#-----PYMOL INSTRUCTIONS-------

# 1. SET WORKING DIR IN PYMOL
from pymol import stored
import re, sys
import pymol 

# change working dir

# view colors
iterate all, print(color)
iterate all, print(b)

# 2. clear the b factor first
>alter ___, b = 0.0
# Alter: modified 1443 atoms.

# 3. input is "byResidue_ID.txt" file from pdb .txt version
>run data2bfactor.py
>data2b_res ____

# 4. sat. refers to min color  and value refers to max color
>run color_b.py
>color_b ____


--------IF ABOVE METHOD FAILS----------
fetch ____
or dele ____
then fetch ____

# 7. file with b factor repleacement for alpha carbon
>infil = open("CADD_CASP8_4jj7_v2.txt", 'r')
>stored.newb = []
>for line in infil.readlines(): stored.newb.append(float(line))
>infil.close()
>print(stored.newb)

>alter ___, b = 0.0
# 8. alter all atoms with b fac values
#>alter ___, b = stored.newb.pop(0)
>alter ___  and n. CA, b=stored.newb.pop(0) # can use alpha2all.py after this

# SAVE updated b factors of CA with new propt. CADD avg

# color based on new b factors of the alpha carbons
>spectrum b, white gray black, ____

# 10. otherwise put "" around object name
>cmd.spectrum("b", "rainbow", "kras", minimum=13.83, maximum=33.25)
# if you dont enter max min will default auto set Spectrum: range ( 0.00000 to 33.25000)

# IF SPECTRUM FAILS TRY....
# script that extends to all atoms so color matches
>run alpha2all.py
>alphaToAll 4jj7
#>alphaToAll "4jj7"

-------CASP8 again-------
# selection chain A and chain B, right click actions, copy to object, new
>run loadBfactors.py  # edit script to have chain specifc scores
>loadBfacts obj01
# same for chain B
>run alpha2all.py
>alphaToAll obj01
>alphaToAll obj02
# set color
>spectrum b, white gray black, obj01
>spectrum b, white gray black, obj02

# fix cartoon after script alters settings
>cartoon automatic, ____


----PYMOL movie----
MOUSE - set 1 button viewing mode

# store a scene, a view that looks nice
scene 001, store
# go to another zoomed in "orient" view of labeled Cys
scene 002, store

# 5 angstrom on this surrounding inhibitor for active site selection
select active, byres all within 5 of inhibitors

# select all cys in a selected object, copy to new object!
sel cysteines, resn cys 

# preset displacts in action,
# CONSURF ran consurf_new.py on .pdb loaded from results output folder
# Set the color to b factor or used alpha to all then did spectrum

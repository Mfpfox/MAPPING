"""
KRAS and CADD scores figure
website script repo: http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
4/17/19
mfpfox

notes: 
- this uses byResidue.txt to set b factor values
- setting residue values, output on pymol shows "keeping..." for the b 
-factor values corresponding to atoms
- run data2bfactor.py
- data2b_res ____ (5ocg)

Please read below for instructions

contains the functions
   data2b_atom(mol='',data_file='')
   data2b_res(mol='',data_file='')
   data2q_atom(mol='',data_file='')
   data2q_res(mol='',data_file='')
"""

import sys, re
from pymol import stored

comment = re.compile('^\s*$|^\s*#')

def atom_data_extract(data_file):
    """
    Read the specified 'by-atom' data file and extract the data from it
    and store it in parallel dictionaries specifying the data
    and residue names (both with keys of chain and residue number and atom name).
    The data file can contain comment lines starting with "#" (on lines by themselves).
    These comment lines are ignored.
    """
    bdat = {}
    chain = ''

    data_lines = open(data_file, 'r')
    count = 0

    for line in data_lines:
        # ignore comment lines (beginning with a '#') or blank lines
        if not comment.match(line):
            words = line.split()
            count += 1

            # check number of columns of data
            if len(words) >= 5:
                chain = words[0]
                resi = words[1]
                resn = words[2]
                name = words[3]
                if chain == '-':
                    chain = ''
                ID = None
                data = float(words[4])
            elif len(words) == 4:
                resi = words[0]
                resn = words[1]
                name = words[2]
                ID = None
                data = float(words[3])
            elif len(words) == 2:
                # ptraj `atominfluct` output use ID and data
                resi = None
                resn = None
                name = None
                ID = int(words[0])
                data = float(words[1])
            else:
                sys.stderr.write("\nError in reading data files -- check number of columns.\n")
                sys.stderr.write("Number of columns: %d  on line number %d\n" % (len(words),count))
                sys.exit(1)

            if ID is None:
                bdat.setdefault(chain, {}).setdefault(resi, {})[name] = (data, resn)
            else:
                # it seem `resn` is not used in anywhere
                bdat.setdefault(chain, {})[ID] = (data, '')

    return bdat

def residue_data_extract(data_file):
    """
    Read the specified 'by-residue' data file and extract the data from it
    and store it in parallel dictionaries specifying the data
    and residue names (both with keys of chain and residue number).
    The data file can contain comment lines starting with "#" (on lines by themselves).
    These comment lines are ignored.
    """
    bdat = {}
    chain = ''

    data_lines = open(data_file, 'r')
    count = 0

    for line in data_lines:
        # ignore comment lines (beginning with a '#') or blank lines
        if not comment.match(line):
            words = line.split()
            count += 1

            # check number of columns of data
            if len(words) >= 4:
                chain = words[0]
                resi = words[1]
                resn = words[2]
                if chain == '-':
                    chain = ''
                data = float(words[3])
            elif len(words) == 3:
                resi = words[0]
                resn = words[1]
                data = float(words[2])
            elif len(words) == 2:
                resi = words[0]
                data = float(words[1])
                resn = ''
            else:
                sys.stderr.write("\nError in reading data files -- check number of columns.\n")
                sys.stderr.write("Number of columns: %d  on line number %d\n" % (len(words),count))
                sys.exit(1)

            bdat.setdefault(chain, {})[resi] = (data, resn)

    return bdat

###########################################################################################
# for testing purposes:
# if calling this as a program on its own, read the data_file name from
# the command line and run residue_data_extract on it. (does not require
# importing cmd from pymol

if __name__ == '__main__':
    data_file_name = sys.argv[1]
    b_dict = residue_data_extract(data_file_name)
    for chain in sorted(b_dict):
        for resi in sorted(b_dict[chain]):
            b, resn = b_dict[chain][resi]
            print("b-factors %s %s %s %s  new B='%s'" % (pdb_file, chain, resn, resi, b))
    sys.exit()


###########################################################################################
# PyMOL stuff

from pymol import cmd
# data2b_atom A # for the chain
def data2b_atom(mol, data_file='toAtom_cadd.txt', property='b', quiet=0):
    """
DESCRIPTION

    Alters the B-factor data by atom.

USAGE

    data2b_atom <mol>, <data_file>

    where <mol> is the molecular object whose B-factor data you wish to modify
    and <data_file> is a file contain the data (one value for each atom)
    The format of <data_file> should be:

         chain resi resn name data
    or
         resi resn name data

    or
         ID data

    (i.e. "chain" is optional if all atoms are in one chain).
    Lines beginning with '#' are ignored as comments.

    Example data lines:

      A ALA 1 N 3.5
      A ALA 1 CA 3.0
      A ALA 1 CB 3.14159
      A ALA 1 C  6.23
      A ALA 1 O    5.1
      A GLY 2 N 10.3
      A GLY 2 CA 1.714
      A GLY 2 C -0.05
      A GLY 2 O -3.12

SEE ALSO

    data2b_res, data2q_atom, data2q_res
    """

    b_dict = atom_data_extract(data_file)
    quiet = int(quiet) == 1

    def b_lookup(chain, resi, name, ID, b):
        def _lookup(chain, resi, name, ID):
            if resi in b_dict[chain] and isinstance(b_dict[chain][resi], dict):
                return b_dict[chain][resi][name][0]
            else:
                # find data by ID
                return b_dict[chain][int(ID)][0]
        try:
            if not chain in b_dict:
                chain = ''
            b = _lookup(chain, resi, name, ID)
            if not quiet: print('///%s/%s/%s new: %f' % (chain, resi, name, b))
        except KeyError:
            if not quiet: print('///%s/%s/%s keeping: %f' % (chain, resi, name, b))
        return b
    stored.b = b_lookup

    cmd.alter(mol, '%s=stored.b(chain, resi, name, ID, %s)' % (property, property))
    cmd.rebuild()



# hardcoded these options
# change data_file for specific mole
# byResidue_chB_3kjn_casp8.txt
def data2b_res(mol, data_file='byResidue_chA_3kjn_casp8.txt', property='b', quiet=0):

    """
DESCRIPTION

    Alters the B-factor data by residue.

USAGE

    data2b_res <mol>, <data_file>

    where <mol> is the molecular object whose B-factor data you wish to modify
    and <data_file> is a file contain the data (one value for each residue)
    The format of <data_file> should be:

         chain resi resn data
    or
         resi resn data
    or
         resi data

    (i.e. "chain" is optional). Lines beginning with '#' are ignored as comments.

SEE ALSO

    data2b_atom, data2q_atom, data2q_res
    """

    b_dict = residue_data_extract(data_file)
    quiet = int(quiet) == 1

    def b_lookup(chain, resi, name, b):
        try:
            if chain in b_dict:
                b = b_dict[chain][resi][0]
            else:
                b = b_dict[''][resi][0]
            if not quiet: print('///%s/%s/%s new: %f' % (chain, resi, name, b))
        except KeyError:
            if not quiet: print('///%s/%s/%s keeping: %f' % (chain, resi, name, b))
        return b
    stored.b = b_lookup

    cmd.alter(mol, '%s=stored.b(chain, resi, name, %s)' % (property, property))
    cmd.rebuild()

def data2q_atom(mol='',data_file=''):
    """
DESCRIPTION

    Alters the occupancy data by atom.

USAGE

    See data2b_atom

SEE ALSO

    data2q_res, data2b_atom, data2b_res
    """
    data2b_atom(mol, data_file, property='q')

def data2q_res(mol='',data_file=''):
    """
DESCRIPTION

    Alters the occupancy data by residue.

USAGE

    See data2b_res

SEE ALSO

    data2q_atom, data2b_atom, data2b_res
    """
    data2b_res(mol, data_file, property='q')

cmd.extend('data2b_res',data2b_res)
cmd.extend('data2b_atom',data2b_atom)
cmd.extend('data2q_res',data2q_res)
cmd.extend('data2q_atom',data2q_atom)

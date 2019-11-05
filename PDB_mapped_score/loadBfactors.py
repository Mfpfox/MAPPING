from pymol import cmd, stored, math
"""
 From kras pro
 # option 1
#def loadBfacts (mol,startaa=1,source="newBfactors.txt", visual="Y"):
# option 2
def loadBfacts (mol,startaa=1,source="newBfactors_MinimumValueReplace2_13.txt", visual="Y"):
   
    mynotes: position on KRAS that dont match ukb sequence (no CADD score) == 
    1, 12, 151, 153, 165-172
    option 1: change all above positions to 1 value (which makes scale for colors 1-33)
    option 2: change all above positions to miniumn of set values that map to ukb sequence (range 13-33)
    in pymol: 
        change working directory
        run loadBfactors.py
        loadBfacts ____ (5ocg)
    problem: 
        only changing b factor of alpha carbons not side chains

    Replaces B-factors with a list of values contained in a plain txt file

    usage: loadBfacts mol, [startaa, [source, [visual]]]

    mol = any object selection (within one single object though)
    startaa = number of first amino acid in 'new B-factors' file (default=1)
    source = name of the file containing new B-factor values (default=newBfactors.txt)
    visual = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)

    example: loadBfacts 1LVM and chain A

"""
# add file as source
#CADD_CASP8_3KJN_chainA_pos223_374.txt
#CADD_CASP8_3KJN_chainB_pos390_379.txt
#CADD_G6PD_2BH9_pos27_515.txt
def loadBfacts (mol, startaa=223,\
                source="CADD_CASP8_3KJN_chainA_pos223_374.txt", visual="Y"):
    obj=cmd.get_object_list(mol)[0]
    cmd.alter(mol,"b=-1.0") # changing them all to -1 before editing
    inFile = open(source, 'r')
    counter = int(startaa) # change start position to integer
    bfacts = []
    for line in inFile.readlines(): # read each line
        bfact = float(line)
        bfacts.append(bfact) # add b factor read in to the list
        cmd.alter("%s and resi %s and n. CA"%(mol, counter), "b=%s"%bfact) # alter object
        # coloring based on bfactor value of alpha carbon (backbone)
        counter = counter+1

    if visual=="Y":
        #cmd.show_as("cartoon",mol)
        # cmd.cartoon("putty", mol)

        # This setting determines how PyMOL renders the transformation of original values into putty related settings
        # 0 value = # normalized nonlinear scaling
        # 4 value =  normalized linear scaling
        # 7 value = absolute linear scaling from the B factor
        # 7- higher CADD means fatter size of putty
        # 0 and 4 look similar, i dont want thickness to translate info so keeping at 0

        #cmd.set("cartoon_putty_scale_min", min(bfacts),obj)
        #cmd.set("cartoon_putty_scale_max", max(bfacts),obj)
        #cmd.set("cartoon_putty_transform", 0,obj)
        #cmd.set("cartoon_putty_transform", 7,obj)
        #cmd.set("cartoon_putty_transform", 4,obj)

        #cmd.set("cartoon_putty_radius", 0.08,obj)
        # color adjust here: bmr, rainbow, yellow_white_red

        # copied code from https://pymolwiki.org/index.php/Color#B-Factors
        # cmd.spectrum("b", "blue_white_red", "%s and n. CA"%prot, minimum=0, maximum=maxval)
        # cmd.ramp_new("ramp_obj", prot, range=[0, minval, maxval], color="[blue, white, red ]")


        cmd.spectrum("b","white_gray_black", "%s and n. CA"%mol, minimum=min(bfacts), maximum=max(bfacts))

        #cmd.ramp_new("count", obj, [min(bfacts), max(bfacts)], "rainbow")
        cmd.create("ca_obj", "%s and n. CA"%mol)
        #cmd.ramp_new("ramp_obj", mol, range=[min(bfacts), max(bfacts)], color="[blue, white, red]")
        cmd.ramp_new("ramp_obj", obj, range=[min(bfacts), max(bfacts)], color="grayscale")

        #cmd.set("surface_color", "ca_obj", mol) 
        cmd.recolor()
        # scale bar options: afmhot grayscale  object   rainbow  traditional
        #grayable     hot          ocean        sludge     
cmd.extend("loadBfacts", loadBfacts);

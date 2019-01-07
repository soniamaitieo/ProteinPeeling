
import glob
import os
import sys


# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]

import pymol

# Call the function below before using any PyMOL modules.
pymol.finish_launching()

from pymol import cmd

#/!\ ATTENTION LE PROGRAMME SE LANCE DEPUIS Src
#CWD = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
#path = CWD + "/results/" + pdb_code + "/PU" + str(PU + 1) + '/'

#path = "/home/soniamai/Bureau/Protein_Peeling/results/1atn/PU1/"
path= sys.argv[1]
#print (path)

for file in glob.iglob(path + '*.pdb'):
    print(file.split("/")[-1])
    cmd.hide("everything" , "all" )
    cmd.load(file, str((file.split("/")[-1])))
    cmd.show("cartoon" , str((file.split("/")[-1])) )
    cmd.hide("line" , str((file.split("/")[-1])) )
    cmd.zoom(str((file.split("/")[-1])) )
    cmd.center(str((file.split("/")[-1])) )
    cmd.save( path + str((file.split("/")[-1])) + ".png", str((file.split("/")[-1])) , "png")


cmd.show("cartoon" , "all")
cmd.zoom("all")
cmd.center("all")
cmd.save(path + str((file.split("/")[-1].split(".")[0])) + ".png" , "all","png")

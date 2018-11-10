import os,sys

directory = sys.path[0]

for subdir in os.listdir(directory):
    dir = os.path.join(directory, subdir)
    if os.path.isdir( dir ):
        dat_file = os.path.join(subdir, "transferDiags.dat")
        os.system( "cp {} ~/Dropbox/{}.dat".format(dat_file,subdir[3:]) )
        
        dat_file = os.path.join(subdir, "diags.dat")
        os.system( "cp {} ~/Dropbox/totEnergy_{}.dat".format(dat_file,subdir[3:]) )

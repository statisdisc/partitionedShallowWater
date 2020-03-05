import sys
import os
import numpy as np

folder_xyz = os.path.join(sys.path[0], "xyzFiles")
os.system( "rm -rf {}".format(folder_xyz) )
os.makedirs( folder_xyz )

for i in xrange(1001):
    folder_i = os.path.join(sys.path[0], str(i))
    folder_xyz_i = os.path.join(folder_xyz, str(i))
    
    os.makedirs( folder_xyz_i )
    
    os.system( "cp {} {}/".format( os.path.join(folder_i,"*.xyz"), folder_xyz_i ) )
        
# os.system( "cp -r {}/ ~/Dropbox/PhD/2019/04_April/".format(folder_xyz) )
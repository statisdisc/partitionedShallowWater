import numpy as np
import os,sys
    
    
    
execfile(os.path.join(sys.path[0],"run_multiple_resolutions_functions.py"))
execfile( os.path.join( sys.path[0], "dependencies.py" ) )


folder_xyz = os.path.join(sys.path[0], "xyzData")
os.system( "rm -rf {}".format(folder_xyz) )
os.makedirs( folder_xyz )
    
folder_oneCol = os.path.join(folder_xyz, "oneCol")
if not os.path.exists( folder_oneCol ):
    os.makedirs( folder_oneCol )
    
folder_threeCols = os.path.join(folder_xyz, "threeCols")
if not os.path.exists( folder_threeCols ):
    os.makedirs( folder_threeCols )

# folder_initial_profiles = os.path.join(sys.path[0], "profiles")
folder_initial_profiles = "/home/statisdisc/OpenFOAM/statisdisc/run/data/risingBubbleResolved/"

dx = np.array([ 10000 ])
folders = [folder_oneCol, folder_threeCols]
folders = [folder_oneCol]

z_oneColumn = np.arange(-50.,10050.,100.)
z_oneColumn[0] = 0.
z_oneColumn[-1] = 10000.

for i in xrange(len(dx)):
    for j in xrange(len(folders)):
        
        #Grid properties
        if folders[j] == folder_oneCol:
            xmax = max( 5, int(0.5*dx[i]/1000.) )
            xmin = -xmax
        if folders[j] == folder_threeCols:
            xmax = max( 5, int(1.5*dx[i]/1000.) )
            xmin = -xmax
        nx = int( round( 1000*(xmax-xmin)/dx[i] ) )
        
        write_blockMeshDict(xmin, xmax, nx)
        
        for k in xrange(20,420,20):
            id = str(k)
            folder_k = os.path.join(folder_initial_profiles, id)
            make_field_files(id, folder_k, z_oneColumn)
        
            os.system( "./run.sh" )
            
            folder_testCase = os.path.join(folders[j], "dx_{}m".format(dx[i]))
            if not os.path.exists( folder_testCase ):
                os.makedirs( folder_testCase )
                
            folder_testCase = os.path.join(folder_testCase, str(k+1))
            if not os.path.exists( folder_testCase ):
                os.makedirs( folder_testCase )
        
            folder_data = os.path.join(sys.path[0], "1")
            os.system( "cp {} {}/".format( os.path.join(folder_data,"*.xyz"), folder_testCase ) )
        
os.system( "cp -r {}/ ~/Dropbox/PhD/2019/04_April/".format(folder_xyz) )
    




    
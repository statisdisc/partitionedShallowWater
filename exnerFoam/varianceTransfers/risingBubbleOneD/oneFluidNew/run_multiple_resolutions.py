import sys
import os
import numpy as np

folder_xyz = os.path.join(sys.path[0], "xyzData_oneFluid")
os.system( "rm -rf {}".format(folder_xyz) )
os.makedirs( folder_xyz )
    
folder_oneCol = os.path.join(folder_xyz, "oneCol")
if not os.path.exists( folder_oneCol ):
    os.makedirs( folder_oneCol )
    
folder_threeCols = os.path.join(folder_xyz, "threeCols")
if not os.path.exists( folder_threeCols ):
    os.makedirs( folder_threeCols )

execfile(os.path.join(sys.path[0],"run_multiple_resolutions_functions.py"))

dx = np.array([ 100, 200, 400, 1000, 2000, 3333, 10000, 50000, 100000, 200000 ])
dx = np.array([ 100, 200, 400, 1000, 2000, 4000, 6666, 20000, 50000, 100000, 200000 ]) #[::-1]
#dx = np.array([ 2000, 3333, 10000, 50000, 100000, 200000 ])
#dx = np.array([ 100 ])
#dx = np.array([ 1000 ])
folders = [folder_oneCol, folder_threeCols]
folders = [folder_oneCol]



for i in xrange(len(dx)):
    for j in xrange(len(folders)):
        
        #Grid properties
        if folders[j] == folder_oneCol:
            xmax = max( 10, int(10*dx[i]/10000.) )
            xmin = -xmax
        if folders[j] == folder_threeCols:
            xmax = max( 10, int(30*dx[i]/10000.) )
            xmin = -xmax
        nx = int( round( 1000*(xmax-xmin)/dx[i] ) )
        
        x = np.linspace(xmin,xmax,nx+1)
        x = x + 0.5*dx[i]/1000.
        x = x[:-1]

        zmin = 0
        zmax = 10
        nz = 100
        z = np.linspace(zmin,zmax,nz+1)
        dz = z[1]-z[0]
        z = z + 0.5*dz
        z = z[:-1]

        #Bubble properties
        base_temp = 300
        radius = 2.
        xcentre = 0.
        zcentre = 2.

        write_blockMeshDict(xmin, xmax, nx)

        write_theta_field(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp)
        write_thetaBuoyant_field(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp)
        write_thetaStable_field(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp)

        write_thetaVar_field(x, z, dx[i]/1000., xcentre, zcentre, radius)
        write_thetaVarBuoyant_field(x, z, dx[i]/1000., xcentre, zcentre, radius)
        write_thetaVarStable_field(x, z, dx[i]/1000., xcentre, zcentre, radius)

        write_sigmaBuoyant_field(x, z, dx[i]/1000., xcentre, zcentre, radius)
        write_sigmaStable_field(x, z, dx[i]/1000., xcentre, zcentre, radius)

        os.system( "./run.sh" )
        
        #for time in xrange(0,1010,10):
        for time in xrange(1000,1010,10):
            folder_testCase = os.path.join(folders[j], "dx_{}m".format(dx[i]))
            #folder_testCase = os.path.join(folders[j], "dx_{}m_t{}".format(dx[i],time))
            if not os.path.exists( folder_testCase ):
                os.makedirs( folder_testCase )
        
            folder_time = os.path.join(sys.path[0], str(time))
            os.system( "cp {} {}/".format( os.path.join(folder_time,"*.xyz"), folder_testCase ) )
        os.system( "cp {} {}/".format( os.path.join(sys.path[0],"*.dat"), folder_testCase ) )
        
os.system( "cp -r {}/ $DROPBOX/PhD/2020/".format(folder_xyz) )

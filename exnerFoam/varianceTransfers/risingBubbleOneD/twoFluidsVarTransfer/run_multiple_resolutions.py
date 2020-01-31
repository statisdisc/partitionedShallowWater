import sys
import os
import numpy as np

dx = np.array([ 100, 200, 400, 1000, 2000, 4000, 6666, 20000, 50000, 100000, 200000 ])
# dx = np.array([ 2000, 4000, 6666, 20000, 50000, 100000, 200000 ])
dx = np.array([ 20000, 50000, 100000, 200000 ])
# dx = np.array([ 20000, 50000 ])
# dx = np.array([ 6666, 20000 ])
#dx = np.array([ 20000 ])

sigma_array = np.linspace(0.,0.2,21)
sigma_array = np.array([0.4])


for k in xrange(len(sigma_array)):
    # folder_xyz = os.path.join(sys.path[0], "xyz_sigma_{}".format(str(sigma_array[k]).replace(".","_")))
    folder_xyz = os.path.join(sys.path[0], "xyzData")
    os.system( "rm -rf {}".format(folder_xyz) )
    os.makedirs( folder_xyz )
        
    folder_oneCol = os.path.join(folder_xyz, "oneCol")
    if not os.path.exists( folder_oneCol ):
        os.makedirs( folder_oneCol )
        
    folder_threeCols = os.path.join(folder_xyz, "threeCols")
    if not os.path.exists( folder_threeCols ):
        os.makedirs( folder_threeCols )

    execfile(os.path.join(sys.path[0],"run_multiple_resolutions_functions.py"))

    folders = [folder_oneCol, folder_threeCols]
    # folders = [folder_oneCol]



    for i in xrange(len(dx)):
        for j in xrange(len(folders)):
            
            #Grid properties
            if folders[j] == folder_oneCol:
                xmax = max( 10, int(10*dx[i]/20000.) )
                xmin = -xmax
            if folders[j] == folder_threeCols:
                xmax = max( 30, int(30*dx[i]/20000.) )
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
            
            sigma = min(0.9, sigma_array[k] * 10000./dx[i])
            
            write_theta_field(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp)
            write_theta_fields(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp, sigma)

            # write_sigmaBuoyant_field(x, z, dx[i]/1000., xcentre, zcentre, radius, sigma)
            # write_sigmaStable_field(x, z, dx[i]/1000., xcentre, zcentre, radius, sigma)
            
            write_blockMeshDict(xmin, xmax, nx)

            os.system( "./run.sh" )
            
            folder_testCase = os.path.join(folders[j], "dx_{}m".format(dx[i]))
            if not os.path.exists( folder_testCase ):
                os.makedirs( folder_testCase )
            
            folder_1000 = os.path.join(sys.path[0], "1000")
            os.system( "cp {} {}/".format( os.path.join(folder_1000,"*.xyz"), folder_testCase ) )
            
    # os.system( "cp -r {}/ ~/Dropbox/PhD/2019/04_April/".format(folder_xyz) )
    os.system( "cp -r {}/ $DROPBOX/PhD/2020/".format(folder_xyz) )

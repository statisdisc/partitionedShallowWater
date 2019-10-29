import sys
import os
import numpy as np

dx = np.array([ 100, 200, 400, 1000, 2000, 3333, 10000, 50000, 100000, 200000 ])
dx = np.array([ 2000, 4000, 6666, 20000, 50000, 100000, 200000 ])
dx = np.array([ 6666, 20000, 50000, 100000, 200000 ])
# dx = np.array([ 20000, 50000, 100000, 200000 ])
# dx = np.array([ 20000, 50000 ])
# dx = np.array([ 20000 ])

dgamma = 0.005
gamma_array = np.arange(0.005,0.1+dgamma,dgamma)
# x_sigma_lim = 1.2
x_sigma_lim = 2.



#Domain properties
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

testcases = []
testcases.append(
    {
        "folder": "01_{}_uniformSigma_0{}",
        "divTransfer": "false",
        "wZeroTransfer": "false",
        "wVarTransfer": "false",
        "directVarianceTransfer": "false",
        "wVarProduction": "false"
    }
)
testcases.append(
    {
        "folder": "02_{}_uniformSigma_0{}_divu",
        "divTransfer": "true",
        "wZeroTransfer": "false",
        "wVarTransfer": "false",
        "directVarianceTransfer": "false",
        "wVarProduction": "false"
    }
)
testcases.append(
    {
        "folder": "03_{}_uniformSigma_0{}_divu_zeroTransfer",
        "divTransfer": "true",
        "wZeroTransfer": "true",
        "wVarTransfer": "false",
        "directVarianceTransfer": "false",
        "wVarProduction": "false"
    }
)
testcases.append(
    {
        "folder": "04_{}_uniformSigma_0{}_divu_wVarTransfer",
        "divTransfer": "true",
        "wZeroTransfer": "false",
        "wVarTransfer": "true",
        "directVarianceTransfer": "true",
        "wVarProduction": "false"
    }
)
testcases.append(
    {
        "folder": "05_{}_uniformSigma_0{}_divu_wVarTransferProduction",
        "divTransfer": "true",
        "wZeroTransfer": "false",
        "wVarTransfer": "true",
        "directVarianceTransfer": "true",
        "wVarProduction": "true"
    }
)
testcases.append(
    {
        "folder": "06_{}_uniformSigma_0{}_divu_wVarCorr",
        "divTransfer": "true",
        "wZeroTransfer": "false",
        "wVarTransfer": "true",
        "directVarianceTransfer": "false",
        "wVarProduction": "false"
    }
)
testcases.append(
    {
        "folder": "07_{}_uniformSigma_0{}_divu_wVarTransferCorr",
        "divTransfer": "true",
        "wZeroTransfer": "false",
        "wVarTransfer": "true",
        "directVarianceTransfer": "false",
        "wVarProduction": "true"
    }
)

for testcase in testcases:

    for k in xrange(len(gamma_array)):
        gamma = gamma_array[k]
        k_string = str(k+1)
        if len(k_string) == 1:
            k_string = "0" + k_string
        folder_basename = testcase["folder"].format( k_string, str(x_sigma_lim).replace(".","") )
        
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
        folders = [folder_oneCol]



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
                
                x = np.linspace(xmin,xmax,nx+1)
                x = x + 0.5*dx[i]/1000.
                x = x[:-1]
                
                write_theta_field(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp)
                write_theta_fields(x, z, dx[i]/1000., xcentre, zcentre, radius, base_temp, x_sigma_lim)
                
                write_blockMeshDict(xmin, xmax, nx)
                write_transferPropertiesDict(gamma, testcase["divTransfer"], testcase["wZeroTransfer"], testcase["wVarTransfer"], testcase["directVarianceTransfer"], testcase["wVarProduction"])

                os.system( "./run.sh" )
                
                folder_testCase = os.path.join(folders[j], "dx_{}m".format(dx[i]))
                if not os.path.exists( folder_testCase ):
                    os.makedirs( folder_testCase )
                
                folder_1000 = os.path.join(sys.path[0], "1000")
                os.system( "cp {} {}/".format( os.path.join(folder_1000,"*.xyz"), folder_testCase ) )
                
        # os.system( "cp -r {}/ $DROPBOX/PhD/2019/".format(folder_xyz) )
        os.system( "mv {}/ $DROPBOX/PhD/2019/{}_gamma_{}".format(folder_xyz, folder_basename, str(gamma).replace(".","")) )
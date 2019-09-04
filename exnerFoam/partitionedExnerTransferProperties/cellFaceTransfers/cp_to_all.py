import os,sys

folders = sorted([d for d in os.listdir( sys.path[0] ) if os.path.isdir(d)])
for folder in folders:

    if folder != "dataFolder":
        #os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/run_identical.sh {}/".format(folder) )
        os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/run_different.sh {}/".format(folder) )
        os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/system/controlDict {}/system/".format(folder) )
        os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/system/transferProperties {}/system/".format(folder) )
        os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/init_0/sigma.* {}/init_0/".format(folder) )


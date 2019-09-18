import os,sys

folders = sorted([d for d in os.listdir( sys.path[0] ) if os.path.isdir(d)])
for folder in folders:

    #if folder != "dataFolder":
    if folder != "dataFolder" and folder != "00_oneFluid":
        #os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/run_identical.sh {}/".format(folder) )
        #os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/run_different.sh {}/".format(folder) )
        #os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/system/controlDict {}/system/".format(folder) )
        os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/system/fvSchemes {}/system/".format(folder) )
        os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/system/setFieldsDict {}/system/".format(folder) )
        #os.system( "cp 01_Qn_Rn_massExplicit_thetaImplicit/init_0/sigma.* {}/init_0/".format(folder) )
        
        os.system( "cp {}/system/transferProperties_sigmaBubbleTest {}/system/transferProperties".format(folder,folder) )
        #os.system( "cp {}/system/transferProperties_minSigmaTest {}/system/transferProperties".format(folder,folder) )


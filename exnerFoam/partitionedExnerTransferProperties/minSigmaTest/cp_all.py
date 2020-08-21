import os,sys

directory = sys.path[0]

data_folder = os.path.join(directory, "dataFolder")
if not os.path.exists(data_folder):
    os.makedirs(data_folder)
    
folders = sorted([d for d in os.listdir( sys.path[0] ) if os.path.isdir(d)])
for folder in folders:

    if folder != "dataFolder":
        testcase_folder = os.path.join(sys.path[0], folder)
        new_folder = os.path.join(data_folder, folder)
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)
            
        os.system( "cp {}/*.dat {}/".format(testcase_folder, new_folder) )
        os.system( "cp {}/1000/*.xyz {}/".format(testcase_folder, new_folder) )

os.system( "cp -r {}/ $DROPBOX/PhD/2019/".format(data_folder) )

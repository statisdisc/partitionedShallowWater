import sys,os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        
folders = sorted([d for d in os.listdir( sys.path[0] ) if os.path.isdir(d)])

print os.getcwd()
for folder in folders[1:]:
    with cd( os.path.join(sys.path[0], folder) ):
        os.system("./run_different.sh")
        os.system("writeCellDataxyz theta")
        os.system("writeCellDataxyz u")
        #os.system("./run_identical.sh")
        
        #os.system("python cp_all.py")

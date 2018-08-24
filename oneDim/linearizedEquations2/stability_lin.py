import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os,sys

    
#Plot stability regions for a range of parameter space.    
def stability_space(h0,h1,c0,Cd,dt,dx,resolution,filename,folder,annotations=False,color=True):
    
    #Set velocities.
    u1 = np.linspace(-1.1*dx/dt,1.1*dx/dt,2*resolution+1)
    u1 += (u1[1] - u1[0])/2.
    if c0 == -10.:
        u0 = u1.copy()
    else:
        u0 = c0 * dx/dt * np.ones(len(u1))
    # print u1*dt/dx
    #Set y-axis to be sqrt(gh) * Delta t / Delta x
    y = np.linspace(0.,1.1,resolution)
    g = y**2 * dx**2/dt**2 / (h0+h1)
    
    dat_filename = os.path.join(sys.path[0], "data_lin_{}res_{}.dat".format(resolution,filename))
    
    if os.path.isfile( dat_filename ):
        stability = np.loadtxt(dat_filename)
    else:
        sys_folder = os.path.join(folder,"system")
        
        dict1_filename = os.path.join(sys_folder,"setFieldsDict")
        dict1_file = open(dict1_filename+"Init", "r+")
        dict1 = dict1_file.read()
        dict1_file.close()
        
        dict2_filename = os.path.join(sys_folder,"fvSolution")
        dict2_file = open(dict2_filename+"Init", "r+")
        dict2 = dict2_file.read()
        dict2_file.close()
    
        stability = np.zeros((len(g),len(u1)))
        for i in xrange(len(stability)):
            print i
            for j in xrange(len(stability[0])):
                h0 = 0.5
                h1 = 0.5
                u0[j] = 0.005 * dx/dt 
                u1[j] = 0.1 * dx/dt 
                g[i] = 0.1**2 * dx**2/dt**2 / (h0+h1)
                
            
                dict1_file = open(dict1_filename, "wb")
                dict1_file.write(dict1 % (h0,h1,u0[j],u1[j], -0.1*h0,-0.1*h1,-0.1*u0[j],-0.1*u1[j]))
                dict1_file.close()
                
                dict2_file = open(dict2_filename, "wb")
                dict2_file.write(dict2 % (g[i]))
                dict2_file.close()
                
                os.system( os.path.join(folder, "allRun.sh") )
                break
                
                delta_e_filename = os.path.join(folder, "energy.dat")
                delta_e_file = open(delta_e_filename, "r+")
                delta_e = delta_e_file.read()
                delta_e_file.close()
                delta_e = delta_e.split("\n")[-2]
                print delta_e
                t = float(delta_e.split(" ")[0])
                delta_e = float(delta_e.split(" ")[-1])
                print delta_e
                
                if float(delta_e) > 0. or t != 5.:
                    stability[i][j] = 0.
                else:
                    stability[i][j] = 1.
            break
    #Adjust axes to center points visually.
    u1 -= (u1[1] - u1[0])/2.
    x = u1 * dt/dx
    y -= (y[1] - y[0])/2.
    
    
    np.savetxt(dat_filename,stability)
    
    ##########################
    #####     PLOTS     ######
    ##########################
    
    X = np.linspace(-1.,1.,50)
    Y = np.sqrt(1 - 1.5*np.sqrt(h0*h1)) * np.abs(X-c0)
    
    #Both stable.
    ss = (1.,1.,1.)
    #Analytic stable, numerical MAYBE stable.
    se = (1.,1.,0.)
    #Analytic MAYBE stable, numerical stable.
    es = (0.,1.,1.)
    #Analytic stable, numerical unstable
    su = (1.,0.4,0.01)
    #Analytic unstable, numerical stable.
    us = (0.,0.,1.)
    #Both MAYBE stable
    ee = (0.7,0.7,0.7)
    #Analytic MAYBE stable, numerical unstable
    eu = (1.,0.,0.)
    #Analytic unstable, numerical MAYBE stable.
    ue = (0.3,0.,0.6)
    #Both unstable.
    uu = (0.4,0.4,0.4)
    
    #Custom colormap
    if color == True:
        cmap_colors = [ue,eu,ee,us,su,es,se]
        # cmap_colors = [us,su,ss,us,su,ss,ss]
        cm = LinearSegmentedColormap.from_list("just_grey",cmap_colors,7)
        
        # cmap_colors = [(1.,0.,0.),(1.,0.4,0.01),(0.,0.,1.),(0.,1.,1.),(0.,1.,0.)]
        # cmap_colors = [(1.,0.,0.),(1.,0.4,0.01),(0.,0.,1.),(0.,0.,1.),(0.,1.,0.)]
        # cmap_colors = [(1.,0.,0.),(1.,0.4,0.01),'#888888','#888888','#888888']
        # cm = LinearSegmentedColormap.from_list("just_grey",cmap_colors,5)
    elif color == False:
        cmap_colors = [(0.9,0.9,0.9),(1.,0.4,0.01),(0.,0.,1.),(0.,0.,1.),(0.,1.,0.)]
        cm = LinearSegmentedColormap.from_list("just_grey",cmap_colors,5)
    
    x_color = "black"
    y_color = "black"
    for i in xrange(4):
        if i == 1:
            x_color = "white"
        elif i == 2:
            x_color = "black"
            y_color = "white"
        elif i == 3:
            x_color = "white"
        
        plt.figure(figsize=(12,6))
        # cs = plt.contourf(x,y,stability_upwind,levels=np.linspace(0., maximum, 11), cmap=cm, extend="both")
        cs = plt.pcolor(x,y,stability,cmap=cm, vmin=0.15, vmax=0.85)
        #Plot expected distribution for 1 fluid.
        plt.plot(np.linspace(-1.,1.,3.),1.-np.abs(np.linspace(-1.,1.,3.)),"k--",linewidth=2.)
        # plt.plot(X,Y,"k--",linewidth=2.)
        cs.cmap.set_under(uu)
        cs.cmap.set_over(ss)
        # cs.set_clim(0., maximum)
        # cb = plt.colorbar(cs)
        
        if annotations == True:
            labels_x = []
            labels_y = []
            if filename == "c0_005_h0_05":
                labels = ["P1","P2","P3"]
                labels_x = [0.1,0.3,0.5]
                labels_y = [0.1,0.1,0.1]
            if filename == "c0_0505_h0_001":
                labels = ["P4","P5"]
                labels_x = [0.2,0.4]
                labels_y = [0.5,0.5]
            if len(labels_x) > 0:
                plt.plot(labels_x,labels_y,"kx",markersize=14)
                for j in xrange(len(labels)):
                    plt.annotate(labels[j], (labels_x[j] - 0.03,labels_y[j] + 0.03),fontsize=15)
        
        plt.xlabel("$\overline{c}_1 = \overline{u}_1 \\frac{\Delta t}{\Delta x}$",fontsize=24,color=x_color)
        plt.ylabel("$\overline{c}_g = \\sqrt{g(\overline{\\eta}_0  + \overline{\\eta}_1)} \\frac{\Delta t}{\Delta x}$",fontsize=24,color=y_color)
        plt.tick_params(axis="both",labelsize=16)
        plt.xlim(x[0],x[-1])
        plt.ylim(y[0] + (y[1]-y[0])/2.,y[-1])
        if c0 == -10:
            # plt.title("$\\eta_0= %s $m, $\\eta_1 = %s$m, $c_0 = c_1$, $C_D = %s$m$^{-1}$s" % (h0,h1,Cd))
            plt.suptitle("$\overline{\\eta}_0= %s $m, $\overline{\\eta}_1 = %s$m, $\overline{c}_0 = \overline{c}_1$" % (h0,h1),fontsize=24)
        else:
            # plt.title("$\\eta_0= %s $m, $\\eta_1 = %s$m, $c_0 = %s$, $C_D = %s$m$^{-1}$s" % (h0,h1,c0,Cd))
            plt.suptitle("$\overline{\\eta}_0= %s $m, $\overline{\\eta}_1 = %s$m, $\overline{c}_0 = %s$" % (h0,h1,c0),fontsize=24)
        plt.gca().set_aspect("equal")
        plt.savefig( os.path.join( sys.path[0], "stability{}_{}.png".format(i+1,filename) ) , bbox_inches='tight' )
        plt.close()
    
    
def main():
    # custom_legend()
    
    folder = sys.path[0]
    
    dt = 0.01
    dx = 1./40.
    
    resolution = 110
    resolution = 55
    # resolution = 25
    # resolution = 10
    
    #Values of height (0), Courant number and drag coefficient.
    #c0 = -10 is used for the setting for c0 = c1.
    h0 = [0.001,0.1,0.5,0.9,0.999]
    h0 = [0.001,0.5,0.999]
    # h0 = [0.,1.]
    # h0 = [0.999]
    c0 = [-10.,0.005,0.505]
    # c0 = [0.005,0.505]
    # c0 = [-10.]
    Cd = [0.,0.1,1.]
    Cd = [0.]
    
    for c in c0:
        for h in h0:
            for d in Cd:
                print "c0 = {}, h0 = {}, Cd = {}".format(c,h,d)
                filename = "c0_{}_h0_{}".format(c,h)
                filename = filename.replace(".0","").replace("1.","1").replace("0.","0").replace("c0_-10_","")
                # stability_space(h,1.-h,c,d,dt,dx,resolution,filename,folder,annotations=False,color=True)
    
    #Single run.
    h = 0.             #h0
    c = 0.005             #c0
    d = 0.              #Cd
    filename = "c0_{}_h0_{}".format(c,h)
    filename = filename.replace(".0","").replace("1.","1").replace("0.","0").replace("c0_-10_","")
    stability_space(h,1.-h,c,d,dt,dx,resolution,filename,folder,annotations=False,color=True)
    
main()
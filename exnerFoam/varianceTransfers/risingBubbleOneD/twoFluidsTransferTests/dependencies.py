def resort_array(arr):
    arr = np.array(arr)
    arr_new = np.zeros((len(arr[0]),len(arr)))
    for i in xrange(len(arr_new)):
        arr_new[i] = arr[:,i]
    return arr_new
    
def read_xyz_1d(filename,folder=""):
    if folder == "":
        folder = sys.path[0]
        
    filename = os.path.join(folder, filename)
    
    #Read file and split into rows
    file = open(filename, "r+")
    string = file.read()
    file.close()
    
    rows = string.split("\n")
    
    data = []
    
    #Remove comment rows and sort into float data columns.
    for row in rows:
        if row != "":
            if row[0] != "#":
                data_row = row.split("  ")
                for i in xrange(len(data_row)):
                    data_row[i] = float(data_row[i])
                    
                data.append(data_row)
            
    #Flip dimensions of data array so that 1 row of "data" is a single variable
    data_temp = np.array(data)
    data = np.zeros((len(data_temp[0]),len(data_temp)))
    for i in xrange(len(data)):
        data[i] = data_temp[:,i]
        
    #Sort all variables by z-axis position (data[2]).
    index_array = data[2].argsort()
    for i in xrange(len(data)):
        data[i] = data[i][index_array[::-1]]
    
    return data
    
def save_colorbar(my_cmap,cmap_min,cmap_max,label,filename,color_over="",color_under=""):
    a = np.array([[cmap_min,cmap_max]])

    plt.figure(figsize=(12, 0.5))
    img = plt.imshow(a, cmap=my_cmap)
    plt.gca().set_visible(False)
    cax = plt.axes([0.1, 0.5, 0.8, 0.3])
    if color_over == "" and color_under == "":
        cbar = plt.colorbar(orientation="horizontal", cax=cax)
    else:
        cbar = plt.colorbar(orientation="horizontal", cax=cax, extend="both")
    cbar.ax.set_xlabel(label,size=16)
    if color_under != "":
        cbar.cmap.set_under(color_under)
    if color_over != "":
        cbar.cmap.set_over(color_over)
    plt.savefig(os.path.join(sys.path[0],"colorbar_horizontal_{}.png".format(filename)))
    plt.close()

    # plt.figure(figsize=(1, 12))
    plt.figure(figsize=(1, 6))
    img = plt.imshow(a, cmap=my_cmap)
    plt.gca().set_visible(False)
    cax = plt.axes([0.1, 0.025, 0.3, 0.95])
    if color_over == "" and color_under == "":
        cbar = plt.colorbar(orientation="vertical", cax=cax)
    else:
        cbar = plt.colorbar(orientation="vertical", cax=cax, extend="both")
    cbar.ax.set_xlabel(label,size=16)
    if color_under != "":
        cbar.cmap.set_under(color_under)
    if color_over != "":
        cbar.cmap.set_over(color_over)
    plt.savefig(os.path.join(sys.path[0],"colorbar_vertical_{}.png".format(filename)))
    plt.close()
    
def bin_data(z, z_ref, data):
    binned_data = np.zeros(len(z_ref))
    binned_data_tally = np.zeros(len(z_ref))
    
    for k in xrange(len(z)):
        dz = z_ref - z[k]
        z_index = np.argmin( np.abs(dz) )
        
        binned_data[z_index] += data[k]
        binned_data_tally[z_index] += 1
        
    binned_data *= 1./binned_data_tally
    
    return binned_data
    
def bin_data_sigma(z, z_ref, data, sigma):
    binned_data = np.zeros(len(z_ref))
    binned_data_tally = np.zeros(len(z_ref))
    
    for k in xrange(len(z)):
        dz = z_ref - z[k]
        z_index = np.argmin( np.abs(dz) )
        
        binned_data[z_index] += sigma[k]*data[k]
        binned_data_tally[z_index] += sigma[k]
        
    for k in xrange(len(binned_data_tally)):
        binned_data_tally[k] = max( binned_data_tally[k], 1e-16 )
        
    binned_data *= 1./binned_data_tally
    
    return binned_data
    
def get_resolved_variance(z, z_ref, data, sigma):
    data_mean = bin_data_sigma(z, z_ref, data, sigma)
    
    binned_data = np.zeros(len(z_ref))
    binned_data_tally = np.zeros(len(z_ref))
    
    for k in xrange(len(z)):
        dz = z_ref - z[k]
        z_index = np.argmin( np.abs(dz) )
        
        binned_data[z_index] += sigma[k]*(data[k]-data_mean[z_index])**2
        binned_data_tally[z_index] += sigma[k]
        
    for k in xrange(len(binned_data_tally)):
        binned_data_tally[k] = max( binned_data_tally[k], 1e-16 )
        
    binned_data *= 1./binned_data_tally
    
    return binned_data
    
def save_legend(gca,name,ncols=4):
    legend_fig = plt.figure(figsize=(26,26))
    legend = plt.figlegend(*gca.get_legend_handles_labels(), loc='center', ncol=ncols, frameon=False)
    # legend = plt.legend(handles=gca, loc='center', ncol=ncols, frameon=False)
    legend_fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(legend_fig.dpi_scale_trans.inverted())
    ll, ur = bbox.get_points()
    x0, y0 = ll
    x1, y1 = ur
    w, h = x1 - x0, y1 - y0
    x1, y1 = x0 + w * 1.1, y0 + h * 1.1
    bbox = tr.Bbox(np.array(((x0, y0),(x1, y1))))
    legend_fig.savefig(os.path.join(sys.path[0],'legend_%s.png' % (name)),bbox_inches=bbox)
    plt.close()
    
def write_field(name, dimensions, field):

    string1 = '''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  dev                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      '''+name+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      '''+dimensions+''';

internalField   nonuniform List<scalar> 
'''+str(len(field))+'''
(
'''

    string2 = ''')
;

boundaryField
{
    ground
    {
        type            zeroGradient;
    }
    top
    {
        type            zeroGradient;
    }
    left
    {
        type            zeroGradient;
    }
    right
    {
        type            zeroGradient;
    }
    defaultFaces
    {
        type            empty;
    }
}


// ************************************************************************* //'''

    folder = os.path.join( sys.path[0], "init_0" )
    file = open(os.path.join(folder,name),"wb")
    file.write(string1)
    for i in xrange(len(field)):
        file.write(str(field[i]) + "\n")
    file.write(string2)
    file.close()
    
def write_vector_field(name, dimensions, u, v, w):

    string1 = '''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  dev                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      '''+name+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      '''+dimensions+''';

internalField   nonuniform List<vector> 
'''+str(len(u))+'''
(
'''

    string2 = ''')
;

boundaryField
{
    ground
    {
        type            zeroGradient;
    }
    top
    {
        type            zeroGradient;
    }
    left
    {
        type            zeroGradient;
    }
    right
    {
        type            zeroGradient;
    }
    defaultFaces
    {
        type            empty;
    }
}


// ************************************************************************* //'''

    folder = os.path.join( sys.path[0], "init_0" )
    file = open(os.path.join(folder,name),"wb")
    file.write(string1)
    for i in xrange(len(u)):
        vector = "(" + str(u[i]) + " " + str(v[i]) + " " + str(w[i]) + ")"
        file.write(vector + "\n")
    file.write(string2)
    file.close()
    
    
    
    
def write_exner_field(name, field):

    foam_string = '''/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  dev
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "1";
    object      Exner;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
%s
(
%s
)
;

boundaryField
{
    ground
    {
        type            partitionedHydrostaticExner;
        gradient        uniform 3.25373e-05;
        value           uniform %s;
    }
    top
    {
        type            partitionedHydrostaticExner;
        gradient        uniform -3.25373e-05;
        value           uniform %s;
    }
    left
    {
        type            partitionedHydrostaticExner;
        gradient        uniform 0;
        value           nonuniform List<scalar>
%s
(
%s
)
;
    }
    right
    {
        type            partitionedHydrostaticExner;
        gradient        uniform 0;
        value           nonuniform List<scalar>
%s
(
%s
)
;
    }
    defaultFaces
    {
        type            empty;
    }
}


// ************************************************************************* //'''
    field_string = ""
    for i in xrange(len(field)):
        field_string += str(field[i]) + "\n"
    
    foam_string = foam_string % (len(field), field_string, str(field[0]), str(field[-1]), len(field), field_string, len(field), field_string)
    
    folder = os.path.join( sys.path[0], "init_0" )
    file = open(os.path.join(folder,name),"wb")
    file.write( foam_string )
    file.close()
    
    
def make_field_files(id, folder, z_oneColumn):
    x,y,z,u,v,w = read_xyz_1d("u.xyz", folder=folder)
    
    sigma_stable = 1.*np.ones(len(w))
    sigma_buoyant = 0.*np.ones(len(w))
    
    #If fluid is rising, label as buoyant fluid
    w_transition = 0.5
    for k in xrange(len(w)):
        if w[k] > w_transition:
            sigma_buoyant[k] = 1.
            sigma_stable[k] = 0.
        elif w[k] > -w_transition:
            sigma_buoyant[k] = np.sin( (w[k]+w_transition)/(2*w_transition) * np.pi/2. )**2
            sigma_stable[k] = 1. - sigma_buoyant[k]
    
    x,y,z,rho = read_xyz_1d("rho.sigma.stable.xyz", folder=folder)
    x,y,z,theta = read_xyz_1d("theta.xyz", folder=folder)
    x,y,z,exner = read_xyz_1d("Exner.xyz", folder=folder)
    
    sigmaRho_stable = sigma_stable*rho
    sigmaRho_buoyant = sigma_buoyant*rho
    
    # sigma_stable_binned = bin_data(z, z_oneColumn, sigma_stable)
    sigma_buoyant_binned = bin_data(z, z_oneColumn, sigma_buoyant)
    sigma_stable_binned = -sigma_buoyant_binned + 1
    
    exner_binned = bin_data_sigma(z, z_oneColumn, exner, rho)
    
    rho_binned = bin_data(z, z_oneColumn, rho)
    rho_stable_binned = bin_data_sigma(z, z_oneColumn, rho, sigma_stable)
    # rho_stable_binned = bin_data(z, z_oneColumn, rho)
    rho_buoyant_binned = bin_data_sigma(z, z_oneColumn, rho, sigma_buoyant)
    # rho_buoyant_binned = bin_data(z, z_oneColumn, rho)
    
    theta_binned = bin_data_sigma(z, z_oneColumn, theta, rho)
    theta_stable_binned = bin_data_sigma(z, z_oneColumn, theta, sigmaRho_stable)
    theta_buoyant_binned = bin_data_sigma(z, z_oneColumn, theta, sigmaRho_buoyant)
    thetaVar_stable_binned = get_resolved_variance(z, z_oneColumn, theta, sigmaRho_stable)
    thetaVar_buoyant_binned = get_resolved_variance(z, z_oneColumn, theta, sigmaRho_buoyant)
    
    u_binned = np.zeros(len(z_oneColumn))
    v_binned = np.zeros(len(z_oneColumn))
    w_binned = bin_data_sigma(z, z_oneColumn, w, rho)
    u_stable_binned = np.zeros(len(z_oneColumn))
    v_stable_binned = np.zeros(len(z_oneColumn))
    w_stable_binned = bin_data_sigma(z, z_oneColumn, w, sigmaRho_stable)
    u_buoyant_binned = np.zeros(len(z_oneColumn))
    v_buoyant_binned = np.zeros(len(z_oneColumn))
    w_buoyant_binned = bin_data_sigma(z, z_oneColumn, w, sigmaRho_buoyant)
    wVar_stable_binned = get_resolved_variance(z, z_oneColumn, w, sigmaRho_stable)
    wVar_buoyant_binned = get_resolved_variance(z, z_oneColumn, w, sigmaRho_buoyant)
    
    write_exner_field("Exner", exner_binned[1:101])
    
    write_field("sigma.stable", "[0 0 0 0 0 0 0]", sigma_stable_binned[1:101])
    write_field("sigma.buoyant", "[0 0 0 0 0 0 0]", sigma_buoyant_binned[1:101])
    
    write_field("rho.stable", "[1 -3 0 0 0 0 0]", rho_stable_binned[1:101])
    write_field("rho.buoyant", "[1 -3 0 0 0 0 0]", rho_buoyant_binned[1:101])
    
    write_field("theta", "[0 0 0 1 0 0 0]", theta_binned[1:101])
    write_field("theta.stable", "[0 0 0 1 0 0 0]", theta_stable_binned[1:101])
    write_field("theta.buoyant", "[0 0 0 1 0 0 0]", theta_buoyant_binned[1:101])
    write_field("thetaVar.stable", "[0 0 0 2 0 0 0]", thetaVar_stable_binned[1:101])
    write_field("thetaVar.buoyant", "[0 0 0 2 0 0 0]", thetaVar_buoyant_binned[1:101])
    
    write_vector_field("u", "[0 1 -1 0 0 0 0]", u_binned[1:101], v_binned[1:101], w_binned[1:101])
    write_vector_field("u.stable", "[0 1 -1 0 0 0 0]", u_binned[1:101], v_binned[1:101], w_stable_binned[1:101])
    write_vector_field("u.buoyant", "[0 1 -1 0 0 0 0]", u_binned[1:101], v_binned[1:101], w_buoyant_binned[1:101])
    write_field("wVar.stable", "[0 2 -2 0 0 0 0]", wVar_stable_binned[1:101])
    write_field("wVar.buoyant", "[0 2 -2 0 0 0 0]", wVar_buoyant_binned[1:101])
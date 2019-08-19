def bubble2D(x,z,xcentre,zcentre,radius):
    r = np.sqrt( (x-xcentre)**2 + (z-zcentre)**2 )
    if r <= radius:
        L = np.sqrt( ((x-xcentre)/radius)**2 + ((z-zcentre)/radius)**2 )
        return 2*np.cos(0.5*np.pi*L)**2
    else:
        return 0

def bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,sigma,mean="all"):
    if len(x) == 1:
        resolution_x = max( 1, int(100 * (xmax-xmin)/10.) )
    else:
        resolution_x = max( 1, int(100 * (x[-1]-x[0])/(10.*len(x))) )
    resolution_z = 10
    
    dz = z[1] - z[0]
    x = np.linspace(xmin,xmax,resolution_x)
    
    
    
    mean_temp_stable = np.zeros(len(z))
    mean_temp_buoyant = np.zeros(len(z))
    var_temp_stable = np.zeros(len(z))
    var_temp_buoyant = np.zeros(len(z))
    
    for k in xrange(len(z)):
        sigma_index = max(1, int(sigma[k] * resolution_z * len(x)))
        
        temperatures = np.zeros((resolution_z,len(x)))
        for j in xrange(resolution_z):
            
            z_new = z[k] - dz/2. + dz*j*1./resolution_z
            for i in xrange(len(x)):
                value = bubble2D(x[i],z_new,xcentre,zcentre,radius)
                temperatures[j][i] = value
        
        temperatures = temperatures.flatten()
        temperature_indices = np.argsort(temperatures)
        hottest_temps = temperatures[temperature_indices[-sigma_index:]]
        coolest_temps = temperatures[temperature_indices[:-sigma_index]]
        
        mean_temp_stable[k] = np.mean(coolest_temps)
        mean_temp_buoyant[k] = np.mean(hottest_temps)
        var_temp_stable[k] = np.var(coolest_temps)
        var_temp_buoyant[k] = np.var(hottest_temps)
            
    return mean_temp_stable,mean_temp_buoyant,var_temp_stable,var_temp_buoyant
    
def bubble2D_mean(x,z,dx,xcentre,zcentre,radius,sigma,mean="all"):
    mean_stable = np.zeros((len(z),len(x)))
    mean_buoyant = np.zeros((len(z),len(x)))
    var_stable = np.zeros((len(z),len(x)))
    var_buoyant = np.zeros((len(z),len(x)))
    
    for i in xrange(len(x)):
        xmin = x[i] - 0.5*dx
        xmax = x[i] + 0.5*dx
        mean_temp_stable,mean_temp_buoyant,var_temp_stable,var_temp_buoyant = bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,sigma,mean=mean)
        mean_stable[:,i] = mean_temp_stable
        mean_buoyant[:,i] = mean_temp_buoyant
        var_stable[:,i] = var_temp_stable
        var_buoyant[:,i] = var_temp_buoyant
        
    return mean_stable,mean_buoyant,var_stable,var_buoyant

def bubble1D_mean_singleFluid(x,z,xmin,xmax,xcentre,zcentre,radius,mean="all"):
    if len(x) == 1:
        resolution_x = max( 1, int(100 * (xmax-xmin)/10.) )
    else:
        resolution_x = max( 1, int(100 * (x[-1]-x[0])/(10.*len(x))) )
    resolution_z = 10
    
    dz = z[1] - z[0]
    x = np.linspace(xmin,xmax,resolution_x)
    
    bubble = np.zeros(len(z))
    
    for k in xrange(len(z)):
        bubble_mean = 0.
        tally = 0
        
        for j in xrange(resolution_z):
            z_new = z[k] - dz/2. + dz*j*1./resolution_z
            for i in xrange(len(x)):
                value = bubble2D(x[i],z_new,xcentre,zcentre,radius)
                if value != 0. or mean == "all":
                    tally += 1
                    bubble_mean += value
            
        bubble[k] += bubble_mean * 1./max(tally,1)
            
    return bubble
    
def bubble2D_mean_singleFluid(x,z,dx,xcentre,zcentre,radius,mean="all"):
    bubble_mean = np.zeros((len(z),len(x)))
    
    for i in xrange(len(x)):
        xmin = x[i] - 0.5*dx
        xmax = x[i] + 0.5*dx
        bubble_mean[:,i] = bubble1D_mean_singleFluid(x,z,xmin,xmax,xcentre,zcentre,radius,mean=mean)
        
    return bubble_mean
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
'''+str(len(field)*len(field[0]))+'''
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
        for j in xrange(len(field[0])):
            file.write(str(field[i][j]) + "\n")
    file.write(string2)
    file.close()
    
def write_blockMeshDict(xmin, xmax, nx):
    string = '''/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.4                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1000;

xmax  %s;
xmin %s;
nx %s;

vertices
(
    ($xmin 1  0)
    ($xmax 1  0)
    ($xmax 1 10)
    ($xmin 1 10)
    ($xmin 0  0)
    ($xmax 0  0)
    ($xmax 0 10)
    ($xmin 0 10)
);


blocks
(
//    hex (0 1 2 3 4 5 6 7) (50 25 1) simpleGrading (1 1 1)
//    hex (0 1 2 3 4 5 6 7) (100 100 1) simpleGrading (1 1 1)
    hex (0 1 2 3 4 5 6 7) ($nx 100 1) simpleGrading (1 1 1)
);

edges
(
);

patches
(
    wall ground
    (
        (1 5 4 0)
    )
    wall top
    (
        (3 7 6 2)
    )
    wall left
    (
        (0 4 7 3)
    )
	wall right
    (
        (2 6 5 1)
    )
	// empty frontAndBack
    // (
    //     (0 3 2 1)
    //     (4 5 6 7)
    // )
);

mergePatchPairs
(
);

// ************************************************************************* //'''

    folder = os.path.join( sys.path[0], "system" )
    file = open(os.path.join(folder,"blockMeshDict"),"wb")
    file.write(string % (xmax,xmin,nx))
    file.close()
    
def write_sigmaBuoyant_field(x, z, dx, xcentre, zcentre, radius, sigma):
    name = "sigma.buoyant"
    dimensions = "[0 0 0 0 0 0 0]"
    
    field = np.ones((len(z),len(x)))
    for i in xrange(len(x)):
        field[:,i] = sigma
    
    write_field(name, dimensions, field)
    
def write_sigmaStable_field(x, z, dx, xcentre, zcentre, radius, sigma):
    name = "sigma.stable"
    dimensions = "[0 0 0 0 0 0 0]"
    
    field = np.ones((len(z),len(x)))
    for i in xrange(len(x)):
        field[:,i] = 1.-sigma
    
    print field
    
    write_field(name, dimensions, field)

def write_theta_field(x, z, dx, xcentre, zcentre, radius, base_temp):
    name = "theta"
    dimensions = "[0 0 0 1 0 0 0]"
    
    field = base_temp + bubble2D_mean_singleFluid(x,z,dx,xcentre,zcentre,radius,mean="all")
    
    write_field(name, dimensions, field)
    
def write_theta_fields(x, z, dx, xcentre, zcentre, radius, base_temp, sigma):
    name = "theta.buoyant"
    dimensions = "[0 0 0 1 0 0 0]"
    
    mean_stable,mean_buoyant,var_stable,var_buoyant = bubble2D_mean(x,z,dx,xcentre,zcentre,radius,sigma,mean="split")
    mean_stable += base_temp
    mean_buoyant += base_temp
    
    field = base_temp + mean_buoyant
    
    write_field("theta.stable", "[0 0 0 1 0 0 0]", mean_stable)
    write_field("theta.buoyant", "[0 0 0 1 0 0 0]", mean_buoyant)
    write_field("thetaVar.stable", "[0 0 0 2 0 0 0]", var_stable)
    write_field("thetaVar.buoyant", "[0 0 0 2 0 0 0]", var_buoyant)
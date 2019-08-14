def bubble2D(x,z,xcentre,zcentre,radius):
    r = np.sqrt( (x-xcentre)**2 + (z-zcentre)**2 )
    if r <= radius:
        L = np.sqrt( ((x-xcentre)/radius)**2 + ((z-zcentre)/radius)**2 )
        return 2*np.cos(0.5*np.pi*L)**2
    else:
        return 0

def bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,mean="all"):
    if len(x) == 1:
        resolution_x = max( 1, int(200 * (xmax-xmin)/10.) )
    else:
        resolution_x = max( 1, int(200 * (x[-1]-x[0])/(10.*len(x))) )
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
                if (value != 0. and abs(x[i]) <= 1.6) or mean == "all":
                    tally += 1
                    bubble_mean += value
            
        bubble[k] += bubble_mean * 1./max(tally,1)
            
    return bubble
    
def bubble2D_mean(x,z,dx,xcentre,zcentre,radius,mean="all"):
    bubble_mean = np.zeros((len(z),len(x)))
    
    for i in xrange(len(x)):
        xmin = x[i] - 0.5*dx
        xmax = x[i] + 0.5*dx
        bubble_mean[:,i] = bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,mean=mean)
        
    return bubble_mean
    
def bubble1D_var(x,z,xmin,xmax,xcentre,zcentre,radius,mean="all"):
    bubble_mean = bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,mean=mean)
    
    if len(x) == 1:
        resolution_x = max( 1, int(200 * (xmax-xmin)/10.) )
    else:
        resolution_x = max( 1, int(200 * (x[-1]-x[0])/(10.*len(x))) )
    resolution_z = 10
    
    dz = z[1] - z[0]
    x = np.linspace(xmin,xmax,resolution_x)
    
    var = np.zeros(len(z))
    
    for k in xrange(len(z)):
        bubble_var = 0.
        tally = 0
        
        for j in xrange(resolution_z):
            z_new = z[k] - dz/2. + dz*j*1./resolution_z
            for i in xrange(len(x)):
                value = bubble2D(x[i],z_new,xcentre,zcentre,radius)
                if (value != 0. and abs(x[i]) <= 1.6) or mean == "all":
                    tally += 1
                    bubble_var += (value - bubble_mean[k])**2
            
        var[k] += bubble_var * 1./max(tally,1)
            
    return var
    
def bubble2D_var(x,z,dx,xcentre,zcentre,radius,mean="all"):
    bubble_var = np.zeros((len(z),len(x)))
    
    for i in xrange(len(x)):
        xmin = x[i] - 0.5*dx
        xmax = x[i] + 0.5*dx
        bubble_var[:,i] = bubble1D_var(x,z,xmin,xmax,xcentre,zcentre,radius,mean=mean)
        
    return bubble_var
    
def bubble1D_sigma(x,z,xmin,xmax,xcentre,zcentre,radius):
    if len(x) == 1:
        resolution_x = max( 1, int(200 * (xmax-xmin)/10.) )
    else:
        resolution_x = max( 1, int(200 * (x[-1]-x[0])/(10.*len(x))) )
    resolution_z = 10
    
    dz = z[1] - z[0]
    x = np.linspace(xmin,xmax,resolution_x)
    
    bubble = np.zeros(len(z))
    
    for k in xrange(len(z)):
        bubble_sigma = 0.
        
        for j in xrange(resolution_z):
            z_new = z[k] - dz/2. + dz*j*1./resolution_z
            for i in xrange(len(x)):
                if bubble2D(x[i],z_new,xcentre,zcentre,radius) != 0 and abs(x[i]) <= 1.6:
                    bubble_sigma += 1.
            
        bubble[k] = bubble_sigma * 1./(len(x)*resolution_z)
        
        # bubble[k] = max( bubble[k], 0.05)
            
    return bubble
    
def bubble2D_sigma(x,z,dx,xcentre,zcentre,radius):
    bubble_var = np.zeros((len(z),len(x)))
    
    for i in xrange(len(x)):
        xmin = x[i] - 0.5*dx
        xmax = x[i] + 0.5*dx
        bubble_var[:,i] = bubble1D_sigma(x,z,xmin,xmax,xcentre,zcentre,radius)
        
    return bubble_var
    
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
    
def write_sigmaBuoyant_field(x, z, dx, xcentre, zcentre, radius):
    name = "sigma.buoyant"
    dimensions = "[0 0 0 0 0 0 0]"
    
    field = bubble2D_sigma(x,z,dx,xcentre,zcentre,radius)
    field = field.clip(min=1e-9,max=1.-1e-9)
    
    write_field(name, dimensions, field)
    
def write_sigmaStable_field(x, z, dx, xcentre, zcentre, radius):
    name = "sigma.stable"
    dimensions = "[0 0 0 0 0 0 0]"
    
    field = -bubble2D_sigma(x,z,dx,xcentre,zcentre,radius) + 1
    field = field.clip(min=1e-9,max=1.-1e-9)
    
    write_field(name, dimensions, field)
        
def write_theta_field(x, z, dx, xcentre, zcentre, radius, base_temp):
    name = "theta"
    dimensions = "[0 0 0 1 0 0 0]"
    
    field = base_temp + bubble2D_mean(x,z,dx,xcentre,zcentre,radius,mean="all")
    
    write_field(name, dimensions, field)

def write_thetaBuoyant_field(x, z, dx, xcentre, zcentre, radius, base_temp):
    name = "theta.buoyant"
    dimensions = "[0 0 0 1 0 0 0]"
    
    field = base_temp + bubble2D_mean(x,z,dx,xcentre,zcentre,radius,mean="split")
    
    write_field(name, dimensions, field)
    
def write_thetaStable_field(x, z, dx, xcentre, zcentre, radius, base_temp):
    name = "theta.stable"
    dimensions = "[0 0 0 1 0 0 0]"
    
    field = base_temp * np.ones((len(z),len(x)))
    
    write_field(name, dimensions, field)
        
def write_thetaVar_field(x, z, dx, xcentre, zcentre, radius):
    name = "thetaVar"
    dimensions = "[0 0 0 2 0 0 0]"
    
    field = bubble2D_var(x,z,dx,xcentre,zcentre,radius,mean="all")
    
    write_field(name, dimensions, field)

def write_thetaVarBuoyant_field(x, z, dx, xcentre, zcentre, radius):
    name = "thetaVar.buoyant"
    dimensions = "[0 0 0 2 0 0 0]"
    
    field = bubble2D_var(x,z,dx,xcentre,zcentre,radius,mean="split")
    
    write_field(name, dimensions, field)
    
def write_thetaVarStable_field(x, z, dx, xcentre, zcentre, radius):
    name = "thetaVar.stable"
    dimensions = "[0 0 0 2 0 0 0]"
    
    field = np.zeros((len(z),len(x)))
    
    write_field(name, dimensions, field)
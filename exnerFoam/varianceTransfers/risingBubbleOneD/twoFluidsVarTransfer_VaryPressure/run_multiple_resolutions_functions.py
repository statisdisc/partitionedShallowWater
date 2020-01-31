def bubble2D(x,z,xcentre,zcentre,radius):
    r = np.sqrt( (x-xcentre)**2 + (z-zcentre)**2 )
    if r <= radius:
        L = np.sqrt( ((x-xcentre)/radius)**2 + ((z-zcentre)/radius)**2 )
        return 2*np.cos(0.5*np.pi*L)**2
    else:
        return 0

def bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,x_sigma_lim,mean="all"):
    if len(x) == 1:
        resolution_x = max( 1, int(100 * (xmax-xmin)/10.) )
    else:
        resolution_x = max( 1, int(100 * (x[-1]-x[0])/(10.*len(x))) )
    resolution_z = 10
    
    dz = z[1] - z[0]
    x = np.linspace(xmin,xmax,resolution_x)
    
    sigma_buoyant = np.zeros(len(z))
    mean_temp_stable = np.zeros(len(z))
    mean_temp_buoyant = np.zeros(len(z))
    var_temp_stable = np.zeros(len(z))
    var_temp_buoyant = np.zeros(len(z))
    
    for k in xrange(len(z)):
        hottest_temps = []
        coolest_temps = []
        sigma0 = 0
        sigma1 = 0
        
        temperatures = np.zeros((resolution_z,len(x)))
        for j in xrange(resolution_z):
            
            z_new = z[k] - dz/2. + dz*j*1./resolution_z
            for i in xrange(len(x)):
                value = bubble2D(x[i],z_new,xcentre,zcentre,radius)
                # temperatures[j][i] = value
                if np.abs(x[i]) <= x_sigma_lim:
                    hottest_temps.append(value)
                    sigma1 += 1
                else:
                    coolest_temps.append(value)
                    sigma0 += 1

        
        sigma_buoyant[k] = sigma1/float(sigma0+sigma1)
        mean_temp_stable[k] = np.mean( np.array(coolest_temps) )
        mean_temp_buoyant[k] = np.mean( np.array(hottest_temps) )
        var_temp_stable[k] = np.var( np.array(coolest_temps) )
        var_temp_buoyant[k] = np.var( np.array(hottest_temps) )
            
    return sigma_buoyant, mean_temp_stable,mean_temp_buoyant,var_temp_stable,var_temp_buoyant
    
def bubble2D_mean(x,z,dx,xcentre,zcentre,radius,x_sigma_lim,mean="all"):
    sigma1 = np.zeros((len(z),len(x)))
    mean_stable = np.zeros((len(z),len(x)))
    mean_buoyant = np.zeros((len(z),len(x)))
    var_stable = np.zeros((len(z),len(x)))
    var_buoyant = np.zeros((len(z),len(x)))
    
    for i in xrange(len(x)):
        xmin = x[i] - 0.5*dx
        xmax = x[i] + 0.5*dx
        sigma_buoyant, mean_temp_stable,mean_temp_buoyant,var_temp_stable,var_temp_buoyant = bubble1D_mean(x,z,xmin,xmax,xcentre,zcentre,radius,x_sigma_lim,mean=mean)
        sigma1[:,i] = sigma_buoyant
        mean_stable[:,i] = mean_temp_stable
        mean_buoyant[:,i] = mean_temp_buoyant
        var_stable[:,i] = var_temp_stable
        var_buoyant[:,i] = var_temp_buoyant
        
    return sigma1,mean_stable,mean_buoyant,var_stable,var_buoyant

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
    

    

def write_theta_field(x, z, dx, xcentre, zcentre, radius, base_temp):
    name = "theta"
    dimensions = "[0 0 0 1 0 0 0]"
    
    field = base_temp + bubble2D_mean_singleFluid(x,z,dx,xcentre,zcentre,radius,mean="all")
    
    write_field(name, dimensions, field)
    
def write_theta_fields(x, z, dx, xcentre, zcentre, radius, base_temp, x_sigma_lim):
    name = "theta.buoyant"
    dimensions = "[0 0 0 1 0 0 0]"
    
    sigma_buoyant, mean_stable,mean_buoyant,var_stable,var_buoyant = bubble2D_mean(x,z,dx,xcentre,zcentre,radius,x_sigma_lim,mean="split")
    mean_stable += base_temp
    mean_buoyant += base_temp
    
    field = base_temp + mean_buoyant
    
    write_field("sigma.stable", "[0 0 0 0 0 0 0]", -sigma_buoyant + 1.)
    write_field("sigma.buoyant", "[0 0 0 0 0 0 0]", sigma_buoyant)
    write_field("theta.stable", "[0 0 0 1 0 0 0]", mean_stable)
    write_field("theta.buoyant", "[0 0 0 1 0 0 0]", mean_buoyant)
    write_field("thetaVar.stable", "[0 0 0 2 0 0 0]", var_stable)
    write_field("thetaVar.buoyant", "[0 0 0 2 0 0 0]", var_buoyant)
    
    
    
    
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
    
    
    
    
def write_transferPropertiesDict(gamma, divTransfer, wZeroTransfer, wVarTransfer, directVarianceTransfer, wVarProduction, thetaVarTransfer, thetaVarTransferSharp, thetaVarTransferSmooth):
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
    object          environmentalProperties;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cloudRadiusMax     cloudRadiusMax   [0 1 0 0 0]     4000;
cloudRadiusMin     cloudRadiusMin   [0 1 0 0 0]     4000;
dragCoeff                                           0;
dragTransferCoeff                                   0.5;
divCoeff                                            1;
Ksigma Ksigma                       [0 2 -1 0 0]    0;//5e4;//2000;
minSigma                                            1e-9;
// Transfer between partitions based on horizontal divergence
wTransfer                                               false;
divTransfer                                             '''+divTransfer+''';
dragTransfer                                            false;
// Transfer between partitions based on laplacian(theta)
thetaTransfer                                           false;
thetaTransferDiffusivity thetaTransferDiffusivity  [0 2 -1 0 0] 0;//5e4;

directVarianceTransfer                                  '''+directVarianceTransfer+''';
thetaVarTransfer                                        '''+thetaVarTransfer+''';
thetaVarTransferSharp                                   '''+thetaVarTransferSharp+''';
thetaVarTransferSmooth                                  '''+thetaVarTransferSmooth+''';
thetaVarTimescale   thetaVarTimeScale   [0 0 1 0 0]     0.01;
wVarTransfer                                            '''+wVarTransfer+''';
wVarTimescale       wVarTimeScale       [0 0 1 0 0]     0.01;
wZeroTransfer                                           '''+wZeroTransfer+''';
localThetaVarTransfer                                   false;
varMassTransfer                                         false;

wVarProduction                                          '''+wVarProduction+''';
wVarProductionTimescale  wVarProductionTimescale [0 -1 2 0 0]  5;
wVarProductionSeparation                                1.;

//Kw Kw [0 2 -1 0 0] 6e3;//4e3;
Ktheta Ktheta [0 2 -1 0 0] 0;
Kw Kw [0 0 1 0 0] %s;

wVarDiffusion                                           false;
KwVariance          KwVariance          [0 2 -1 0 0]    2e6;

// ************************************************************************* //'''

    folder = os.path.join( sys.path[0], "system" )
    file = open(os.path.join(folder,"transferProperties"),"wb")
    file.write(string % (gamma))
    file.close()

#!/bin/bash -e

#Remove previous simulations.
rm -rf 0 [0-9]*[0-9] constant/polyMesh core log legends gmt.history
#rm *.dat

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
cp init_0/u 0/u.stable
cp init_0/u 0/u.buoyant
cp init_0/Uf 0/Uf.stable
cp init_0/Uf 0/Uf.buoyant
cp init_0/volFlux 0/volFlux.stable
cp init_0/volFlux 0/volFlux.buoyant

#Create initial partition fraction distributions.
cp init_0/sigma.stable 0/
cp init_0/sigma.buoyant 0/
cp init_0/sigmah.stable 0/
cp init_0/sigmah.buoyant 0/
cp init_0/h 0/

cp init_0/U 0/U.stable
cp init_0/U 0/U.buoyant
cp init_0/UF 0/UF.stable
cp init_0/UF 0/UF.buoyant
cp init_0/VOLFLUX 0/VOLFLUX.stable
cp init_0/VOLFLUX 0/VOLFLUX.buoyant
cp init_0/SIGMAH 0/SIGMAH.stable
cp init_0/SIGMAH 0/SIGMAH.buoyant

setFields
#sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

# Plot the initial conditions
export time=0
writeCellDataxyz -time $time buoyant.sigma
gmtPlot plots/plotSigma.gmt

# Solve the SWE
#partitionedShallowWaterFoamFluxExp >& log & sleep 0.01; tail -f log
partitionedLinearizedShallowWaterFoam

# Plots of results
for ((itime=0;itime<51;itime=itime+1))
{
export time=$(bc <<<"scale=1; $itime/10" )
writeCellDataxyz -time $time sigmah.stable
writeCellDataxyz -time $time sigmah.buoyant
writeCellDataxyz -time $time u.stable
writeCellDataxyz -time $time u.buoyant
}

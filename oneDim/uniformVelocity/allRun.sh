#!/bin/bash -e

#Remove previous simulations.
rm -rf 0 [0-9]*[0-9] constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
cp init_0/U 0/stable.U
cp init_0/U 0/buoyant.U
cp init_0/Uf 0/stable.Uf
cp init_0/Uf 0/buoyant.Uf

#Create initial partition fraction distributions.
cp init_0/buoyant.sigma 0/
cp init_0/stable.sigma 0/
setFields
sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

cp init_0/h 0/

# Plot the initial conditions
export time=0
writeCellDataxyz -time $time buoyant.sigma
gmtPlot plots/plotSigma.gmt

# Solve the SWE
#partitionedShallowWaterFoamFluxExp >& log & sleep 0.01; tail -f log
partitionedShallowWaterFoamFluxExp

# Plots of results
#for ((itime=0;itime<3;itime=itime+1))
for ((itime=0;itime<21;itime=itime+2))
{
export time=$(bc <<<"scale=1; $itime/10" )
writeCellDataxyz -time $time buoyant.sigma
writeCellDataxyz -time $time h
writeCellDataxyz -time $time buoyant.Uf
writeCellDataxyz -time $time stable.Uf
writeCellDataxyz -time $time stable.u
# gmtPlot plots/plotSigma.gmt
# gmtPlot plots/plotU.gmt
}

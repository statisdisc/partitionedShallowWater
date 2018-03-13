#!/bin/bash -e

#Remove previous simulations.
rm -rf 0 [0-9]* [0-9]*[0-9] constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
cp init_0/stable.u 0/stable.u
cp init_0/buoyant.u 0/buoyant.u
cp init_0/stable.Uf 0/stable.Uf
cp init_0/buoyant.Uf 0/buoyant.Uf

#Create initial partition fraction distributions.
cp init_0/buoyant.sigma 0/
cp init_0/stable.sigma 0/
sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

cp init_0/h 0/

# Solve the SWE
#partitionedShallowWaterFoamFluxExp >& log & sleep 0.01; tail -f log
partitionedShallowWaterFoamFluxExp

# Plots of results
for ((itime=0;itime<101;itime=itime+1))
{
export time=$(bc <<<"scale=2; $itime/100" )
writeCellDataxyz -time $time buoyant.sigma
writeCellDataxyz -time $time stable.sigma
writeCellDataxyz -time $time h
writeCellDataxyz -time $time buoyant.Uf
writeCellDataxyz -time $time stable.Uf
writeCellDataxyz -time $time stable.flux
writeCellDataxyz -time $time buoyant.flux
}

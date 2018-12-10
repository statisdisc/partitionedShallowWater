#!/bin/bash -e

#Remove previous simulations.
rm -rf [0-9]* [0-9]*[0-9] constant/polyMesh core log legends gmt.history
#rm *.dat

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
cp init/u 0/u.stable
cp init/u 0/u.buoyant
cp init/Uf 0/Uf.stable
cp init/Uf 0/Uf.buoyant
cp init/volFlux 0/volFlux.stable
cp init/volFlux 0/volFlux.buoyant

#Create initial partition fraction distributions.
cp init/sigma.stable 0/
cp init/sigma.buoyant 0/
cp init/h 0/

setFields
#sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

# Plot the initial conditions
export time=0
writeCellDataxyz -time $time buoyant.sigma
gmtPlot plots/plotSigma.gmt

# Solve the SWE
#partitionedShallowWaterFoamFluxExp >& log & sleep 0.01; tail -f log
partitionedShallowWaterFoamAdvExp

# Plots of results
#for ((itime=0;itime<501;itime=itime+1))
#{
#export time=$(bc <<<"scale=10; $itime/100" )
#writeCellDataxyz -time $time sigmah.stable
#writeCellDataxyz -time $time sigmah.buoyant
#writeCellDataxyz -time $time u.stable
#writeCellDataxyz -time $time u.buoyant
#writeCellDataxyz -time $time energy
#}
#find [-0-9]*/ ! -name '*.xyz' -type f -exec rm -f {} +

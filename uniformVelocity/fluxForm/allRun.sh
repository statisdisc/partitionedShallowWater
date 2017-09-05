#!/bin/bash -e

#Remove previous simulations.
rm -rf 0 [0-9]*[0-9] constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
#setVelocityField
#mv 0/U 0/stable.U
#mv 0/Uf 0/stable.Uf
#setVelocityField
#mv 0/U 0/buoyant.U
#mv 0/Uf 0/buoyant.Uf
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

# Solve the SWE
#partitionedShallowWaterFoamAdvExp >& log & sleep 0.01; tail -f log


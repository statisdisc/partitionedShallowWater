#!/bin/bash -e

#Remove previous simulations.
rm -rf 0 [0-9]*[0-9] constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
#setVelocityField
cp init_0/U 0/
cp init_0/Uf 0/

cp init_0/h 0/

# Solve the SWE
#shallowWaterFoam_hu >& log & sleep 0.01; tail -f log

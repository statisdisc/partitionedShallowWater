#!/bin/bash -e

# Create the case, run and post-process

# clear out old stuff
rm -rf [0-9]*  constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir -p 0
cp init_0/* 0
setVelocityField
mv 0/U 0/u

setBalancedHeight
mv 0/u 0/U
mv 0/phi 0/volFlux

# Solve the SWE
shallowWaterFoamAdvExp >& log & sleep 0.01; tail -f log

# Plots
time=10000
for plot in hu sigma; do
    gmtFoam -time $time $plot
    gv $time/$plot.pdf &
done

gmtPlot plots/plotEnergy.gmt

# Plot sigma and hu for all times
for plot in hu sigma; do
    gmtFoam $plot
    eps2gif $plot.gif 0/$plot.pdf ?????/$plot.pdf ??????/$plot.pdf
done


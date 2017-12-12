#!/bin/bash -e

# Create the case, run and post-process

# clear out old stuff
rm -rf 0 [0-9]*[0-9] constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh

#Create velocity field for each partition.
mkdir 0
setVelocityField
cp 0/U 0/stable.u
cp 0/Uf 0/stable.Uf
mv 0/U 0/buoyant.u
mv 0/Uf 0/buoyant.Uf
#cp init_0/U 0/stable.u
#cp init_0/Uf 0/stable.Uf

#Create initial partition fraction distributions.
cp init_0/buoyant.sigma 0/
cp init_0/stable.sigma 0/
#setFields
#sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

#Set permanent transfer term regions between partitions.
cp init_0/stable.source 0/
cp init_0/buoyant.source 0/

cp init_0/h 0/

# Solve the SWE
partitionedShallowWaterFoamAdvExp >& log & sleep 0.01; tail -f log

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


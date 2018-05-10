#!/bin/bash -e

# Create the case, run and post-process

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log legends gmt.history

#Create mesh.
blockMesh


mkdir 0

#Create velocity field for each partition.
cp init_0/U 0/U
cp init_0/Uf 0/Uf
#setVelocityField

#mv 0/U 0/u
cp init_0/h 0/
#setBalancedHeight
#mv 0/u 0/U

cp 0/U 0/stable.u
cp 0/Uf 0/stable.Uf
#cp 0/phi 0/stable.volFlux
mv 0/U 0/buoyant.u
mv 0/Uf 0/buoyant.Uf
#mv 0/phi 0/buoyant.volFlux

cp init_0/stable.volFlux 0/
cp init_0/buoyant.volFlux 0/

#cp init_0/U 0/stable.u
#cp init_0/Uf 0/stable.Uf

#Create initial partition fraction distributions.
cp init_0/buoyant.sigma 0/
cp init_0/stable.sigma 0/
#setFields
#sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

#Set permanent transfer term regions between partitions.
cp init_0/stable.sink 0/
setFields -dict system/stableSinkDict
cp init_0/buoyant.sink 0/
setFields -dict system/buoyantSinkDict

#Set permanent momentum source regions.
cp init_0/stable.momentumSource 0/
cp init_0/buoyant.momentumSource 0/
setFields -dict system/buoyantMomentumSourceDict

#Set permanent momentum decay regions.
cp init_0/stable.momentumSink 0/
cp init_0/buoyant.momentumSink 0/
#setFields -dict system/buoyantMomentumSinkDict

cp init_0/gravity 0/

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


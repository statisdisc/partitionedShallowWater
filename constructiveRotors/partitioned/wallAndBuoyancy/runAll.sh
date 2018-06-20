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

cp 0/U 0/u.stable
cp 0/Uf 0/Uf.stable
mv 0/U 0/u.buoyant
mv 0/Uf 0/Uf.buoyant

cp init_0/volFlux.stable 0/
cp init_0/volFlux.buoyant 0/

#Create initial partition fraction distributions.
cp init_0/sigma.buoyant 0/
cp init_0/sigma.stable 0/

#Set permanent transfer term regions between partitions.
cp init_0/sink.stable 0/
setFields -dict system/stableSinkDict
cp init_0/sink.buoyant 0/

#Set permanent momentum source regions.
cp init_0/momentumSource.stable 0/
cp init_0/momentumSource.buoyant 0/
setFields -dict system/buoyantMomentumSourceDict

#Set permanent momentum decay regions.
cp init_0/momentumSink.stable 0/
cp init_0/momentumSink.buoyant 0/
#setFields -dict system/buoyantMomentumSinkDict

cp init_0/gravity 0/

# Solve the SWE
partitionedShallowWaterFoamAdvExp >& log & sleep 0.01; tail -f log

# Plots
time=10000
for plot in sigmaU sigma; do
    gmtFoam -time $time $plot
    gv $time/$plot.pdf &
done

gmtPlot plots/plotEnergy.gmt

# Plot sigma and hu for all times
for plot in sigmaU; do
    gmtFoam $plot
    eps2gif $plot.gif 0/$plot.pdf ??/$plot.pdf ???/$plot.pdf ????/$plot.pdf
done


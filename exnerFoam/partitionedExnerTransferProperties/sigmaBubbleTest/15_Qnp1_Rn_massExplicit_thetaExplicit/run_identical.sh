#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh
blockMesh

# hydrostatically balanced initial conditions
rm -rf [0-9]* core
mkdir 0
cp -r init_0/* 0
setExnerBalancedH

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedFluxBuoyantExner/partitionedHydrostaticExner/g' 0/Exner

# Add a warm perturnation
cp 0/theta 0/theta_init
#cp constant/initialProperties0 constant/initialProperties
#makeHotBubble
#mv 0/theta 0/stable.theta
#cp constant/initialProperties1 constant/initialProperties

makeHotBubble
cp 0/theta 0/theta.stable
mv 0/theta 0/theta.buoyant

# Partition into stable and buoyant fluids
#cp 0/buoyant.theta 0/stable.theta
#mv 0/theta_init 0/buoyant.theta
mv 0/Uf 0/Uf.stable
cp 0/Uf.stable 0/Uf.buoyant
rm 0/thetaf

# create initial conditions
#setFields
#sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

# Solve Euler equations
#multiFluidFoamAdditionalTransfers >& log & sleep 0.01; tail -f log
multiFluidFoamAdditionalTransfers


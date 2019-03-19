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
cp init_0/theta 0/theta
#makeHotBubble
mv 0/theta 0/theta.stable
cp 0/theta.stable 0/theta.buoyant

cp init_0/thetaVar 0/thetaVar.stable
cp init_0/thetaVar 0/thetaVar.buoyant

# Partition into stable and buoyant fluids
#cp 0/buoyant.theta 0/stable.theta
#mv 0/theta_init 0/buoyant.theta
mv 0/Uf 0/Uf.stable
cp 0/Uf.stable 0/Uf.buoyant

mv 0/u 0/u.stable
cp 0/u.stable 0/u.buoyant

rm 0/thetaf

# create initial conditions
setFields
#sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

# Plot initial conditions
time=0
gmtFoam sigmaTheta -time $time
#gv $time/sigmaTheta.pdf &

# Solve Euler equations
multiFluidFoam
#multiFluidFoam >& log & sleep 0.01; tail -f log

writeCellDataxyz -time 1000 theta
writeCellDataxyz -time 1000 theta.stable
writeCellDataxyz -time 1000 theta.buoyant
writeCellDataxyz -time 1000 thetaVar.stable
writeCellDataxyz -time 1000 thetaVar.buoyant
writeCellDataxyz -time 1000 u.stable
writeCellDataxyz -time 1000 u.buoyant
writeCellDataxyz -time 1000 sigma.stable
writeCellDataxyz -time 1000 sigma.buoyant



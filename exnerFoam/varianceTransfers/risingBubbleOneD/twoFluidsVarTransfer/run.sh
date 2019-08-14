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
#makeHotBubble
cp init_0/theta.stable 0/
cp init_0/theta.buoyant 0/

cp init_0/thetaVar.stable 0/
cp init_0/thetaVar.buoyant 0/

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
#gmtFoam sigmaTheta -time $time
#gv $time/sigmaTheta.pdf &

# Solve Euler equations
#multiFluidFoam >& log & sleep 0.01; tail -f log
multiFluidFoamVarTransfers

writeCellDataxyz theta
writeCellDataxyz theta.stable
writeCellDataxyz theta.buoyant
writeCellDataxyz thetaVar.stable
writeCellDataxyz thetaVar.buoyant
writeCellDataxyz u
writeCellDataxyz u.stable
writeCellDataxyz u.buoyant
writeCellDataxyz wVar.stable
writeCellDataxyz wVar.buoyant
writeCellDataxyz sigma.stable
writeCellDataxyz sigma.buoyant
#writeCellDataxyz -time 1000 sigma.buoyant

# animate the results
#for field in theta sigma thetaU.stable thetaU.buoyant; do
#    gmtFoam $field
#    cp 0/$field.pdf ~/Dropbox/PhD/2019/02_February/05_RisingBubble_ThetaErrors/00_Temp/$field.init.pdf
#    cp 1000/$field.pdf ~/Dropbox/PhD/2019/02_February/05_RisingBubble_ThetaErrors/00_Temp/
#    eps2gif $field.gif 0/$field.pdf ??/$field.pdf ???/$field.pdf ????/$field.pdf
#done

#cp 1000/theta.xyz ~/Dropbox/PhD/2019/02_February/05_RisingBubble_ThetaErrors/00_Temp/
#cp 1000/sigma.buoyant.xyz ~/Dropbox/PhD/2019/02_February/05_RisingBubble_ThetaErrors/00_Temp/
#cp *.gif ~/Dropbox/PhD/2019/02_February/05_RisingBubble_ThetaErrors/00_Temp/


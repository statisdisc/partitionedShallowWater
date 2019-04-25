#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh
blockMesh

# hydrostatically balanced initial conditions
rm -rf [0-9]* core
mkdir 0
cp -r init_0/* 0
#setExnerBalancedH

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedFluxBuoyantExner/partitionedHydrostaticExner/g' 0/Exner

# Add a warm perturnation
cp 0/theta 0/theta_init

rm 0/thetaf


# Solve Euler equations
#multiFluidFoam >& log & sleep 0.01; tail -f log
multiFluidFoamTransferTests

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
writeCellDataxyz rho.sigma.stable
writeCellDataxyz rho.sigma.buoyant
writeCellDataxyz massTransferDivTransfer.stable.buoyant
writeCellDataxyz massTransferDivTransfer.buoyant.stable
writeCellDataxyz massTransferSigmaDiffusion.stable.buoyant
writeCellDataxyz massTransferSigmaDiffusion.buoyant.stable
writeCellDataxyz massTransferThetaDiffusion.stable.buoyant
writeCellDataxyz massTransferThetaDiffusion.buoyant.stable
writeCellDataxyz massTransferThetaVar.stable.buoyant
writeCellDataxyz massTransferThetaVar.buoyant.stable
writeCellDataxyz massTransferWvar.stable.buoyant
writeCellDataxyz massTransferWvar.buoyant.stable
writeCellDataxyz thetaTransferThetaVar.stable.buoyant
writeCellDataxyz thetaTransferThetaVar.buoyant.stable
writeCellDataxyz thetaVarTransferThetaVar.stable.buoyant
writeCellDataxyz thetaVarTransferThetaVar.buoyant.stable
writeCellDataxyz thetaTransferThetaVarAlternative.stable.buoyant
writeCellDataxyz thetaTransferThetaVarAlternative.buoyant.stable
writeCellDataxyz thetaVarTransferThetaVarAlternative.stable.buoyant
writeCellDataxyz thetaVarTransferThetaVarAlternative.buoyant.stable
writeCellDataxyz wTransferWvar.stable.buoyant
writeCellDataxyz wTransferWvar.buoyant.stable
writeCellDataxyz wVarTransferWvar.stable.buoyant
writeCellDataxyz wVarTransferWvar.buoyant.stable
writeCellDataxyz wTransferWvarAlternative.stable.buoyant
writeCellDataxyz wTransferWvarAlternative.buoyant.stable
writeCellDataxyz wVarTransferWvarAlternative.stable.buoyant
writeCellDataxyz wVarTransferWvarAlternative.buoyant.stable
writeCellDataxyz wTransferGradDiv.stable.stable
writeCellDataxyz wTransferGradDiv.buoyant.buoyant
writeCellDataxyz wTransferDrag.stable.buoyant
writeCellDataxyz wTransferDrag.buoyant.stable

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


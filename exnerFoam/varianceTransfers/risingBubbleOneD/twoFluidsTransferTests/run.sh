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
writeCellDataxyz thetaT.stable
writeCellDataxyz thetaT.buoyant
writeCellDataxyz u
writeCellDataxyz u.stable
writeCellDataxyz u.buoyant
writeCellDataxyz wVar.stable
writeCellDataxyz wVar.buoyant
writeCellDataxyz wT.stable
writeCellDataxyz wT.buoyant
writeCellDataxyz sigma.stable
writeCellDataxyz sigma.buoyant
writeCellDataxyz rho.sigma.stable
writeCellDataxyz rho.sigma.buoyant

writeCellDataxyz massTransferAnalytic.stable.buoyant
writeCellDataxyz massTransferAnalytic.buoyant.stable
writeCellDataxyz massTransferDivTransfer.stable.buoyant
writeCellDataxyz massTransferDivTransfer.buoyant.stable
writeCellDataxyz massTransferSigmaDiffusion.stable.buoyant
writeCellDataxyz massTransferSigmaDiffusion.buoyant.stable
writeCellDataxyz massTransferThetaDiffusion.stable.buoyant
writeCellDataxyz massTransferThetaDiffusion.buoyant.stable
writeCellDataxyz massTransferEntrainment.stable.buoyant
writeCellDataxyz massTransferEntrainment.buoyant.stable
                 
writeCellDataxyz thetaTransferAnalyticMean.stable
writeCellDataxyz thetaTransferAnalyticMean.buoyant
writeCellDataxyz thetaTransferAnalyticZero.stable
writeCellDataxyz thetaTransferAnalyticZero.buoyant
writeCellDataxyz thetaTransferAnalyticVar.stable
writeCellDataxyz thetaTransferAnalyticVar.buoyant
writeCellDataxyz thetaTransferAnalyticMeanVar.stable
writeCellDataxyz thetaTransferAnalyticMeanVar.buoyant
#writeCellDataxyz thetaTransferDivMean.stable
#writeCellDataxyz thetaTransferDivMean.buoyant
#writeCellDataxyz thetaTransferDivVar.stable
#writeCellDataxyz thetaTransferDivVar.buoyant
#writeCellDataxyz thetaTransferDivMeanVar.stable
#writeCellDataxyz thetaTransferDivMeanVar.buoyant
#writeCellDataxyz thetaTransferSigmaDiffMean.stable
#writeCellDataxyz thetaTransferSigmaDiffMean.buoyant
#writeCellDataxyz thetaTransferSigmaDiffVar.stable
#writeCellDataxyz thetaTransferSigmaDiffVar.buoyant
#writeCellDataxyz thetaTransferSigmaDiffMeanVar.stable
#writeCellDataxyz thetaTransferSigmaDiffMeanVar.buoyant
#writeCellDataxyz thetaTransferThetaDiffMean.stable
#writeCellDataxyz thetaTransferThetaDiffMean.buoyant
#writeCellDataxyz thetaTransferThetaDiffVar.stable
#writeCellDataxyz thetaTransferThetaDiffVar.buoyant
#writeCellDataxyz thetaTransferThetaDiffMeanVar.stable
#writeCellDataxyz thetaTransferThetaDiffMeanVar.buoyant
                 
writeCellDataxyz thetaVarTransferAnalyticMean.stable
writeCellDataxyz thetaVarTransferAnalyticMean.buoyant
writeCellDataxyz thetaVarTransferAnalyticZero.stable
writeCellDataxyz thetaVarTransferAnalyticZero.buoyant
writeCellDataxyz thetaVarTransferAnalyticVar.stable
writeCellDataxyz thetaVarTransferAnalyticVar.buoyant
writeCellDataxyz thetaVarTransferAnalyticMeanVar.stable
writeCellDataxyz thetaVarTransferAnalyticMeanVar.buoyant
#writeCellDataxyz thetaVarTransferDivMean.stable
#writeCellDataxyz thetaVarTransferDivMean.buoyant
#writeCellDataxyz thetaVarTransferDivVar.stable
#writeCellDataxyz thetaVarTransferDivVar.buoyant
#writeCellDataxyz thetaVarTransferDivMeanVar.stable
#writeCellDataxyz thetaVarTransferDivMeanVar.buoyant
#writeCellDataxyz thetaVarTransferSigmaDiffMean.stable
#writeCellDataxyz thetaVarTransferSigmaDiffMean.buoyant
#writeCellDataxyz thetaVarTransferSigmaDiffVar.stable
#writeCellDataxyz thetaVarTransferSigmaDiffVar.buoyant
#writeCellDataxyz thetaVarTransferSigmaDiffMeanVar.stable
#writeCellDataxyz thetaVarTransferSigmaDiffMeanVar.buoyant
#writeCellDataxyz thetaVarTransferThetaDiffMean.stable
#writeCellDataxyz thetaVarTransferThetaDiffMean.buoyant
#writeCellDataxyz thetaVarTransferThetaDiffVar.stable
#writeCellDataxyz thetaVarTransferThetaDiffVar.buoyant
#writeCellDataxyz thetaVarTransferThetaDiffMeanVar.stable
#writeCellDataxyz thetaVarTransferThetaDiffMeanVar.buoyant

writeCellDataxyz wTransferAnalyticMean.stable
writeCellDataxyz wTransferAnalyticMean.buoyant
writeCellDataxyz wTransferAnalyticVar.stable
writeCellDataxyz wTransferAnalyticVar.buoyant
writeCellDataxyz wTransferAnalyticMeanVar.stable
writeCellDataxyz wTransferAnalyticMeanVar.buoyant
writeCellDataxyz wTransferAnalyticZero.stable
writeCellDataxyz wTransferAnalyticZero.buoyant
#writeCellDataxyz wTransferDivMean.stable
#writeCellDataxyz wTransferDivMean.buoyant
#writeCellDataxyz wTransferDivVar.stable
#writeCellDataxyz wTransferDivVar.buoyant
#writeCellDataxyz wTransferDivMeanVar.stable
#writeCellDataxyz wTransferDivMeanVar.buoyant
#writeCellDataxyz wTransferDivZero.stable
#writeCellDataxyz wTransferDivZero.buoyant
#writeCellDataxyz wTransferSigmaDiffMean.stable
#writeCellDataxyz wTransferSigmaDiffMean.buoyant
#writeCellDataxyz wTransferSigmaDiffVar.stable
#writeCellDataxyz wTransferSigmaDiffVar.buoyant
#writeCellDataxyz wTransferSigmaDiffMeanVar.stable
#writeCellDataxyz wTransferSigmaDiffMeanVar.buoyant
#writeCellDataxyz wTransferSigmaDiffZero.stable
#writeCellDataxyz wTransferSigmaDiffZero.buoyant
#writeCellDataxyz wTransferThetaDiffMean.stable
#writeCellDataxyz wTransferThetaDiffMean.buoyant
#writeCellDataxyz wTransferThetaDiffVar.stable
#writeCellDataxyz wTransferThetaDiffVar.buoyant
#writeCellDataxyz wTransferThetaDiffMeanVar.stable
#writeCellDataxyz wTransferThetaDiffMeanVar.buoyant
#writeCellDataxyz wTransferThetaDiffZero.stable
#writeCellDataxyz wTransferThetaDiffZero.buoyant

writeCellDataxyz wVarTransferAnalyticMean.stable
writeCellDataxyz wVarTransferAnalyticMean.buoyant
writeCellDataxyz wVarTransferAnalyticVar.stable
writeCellDataxyz wVarTransferAnalyticVar.buoyant
writeCellDataxyz wVarTransferAnalyticMeanVar.stable
writeCellDataxyz wVarTransferAnalyticMeanVar.buoyant
writeCellDataxyz wVarTransferAnalyticZero.stable
writeCellDataxyz wVarTransferAnalyticZero.buoyant
#writeCellDataxyz wVarTransferDivMean.stable
#writeCellDataxyz wVarTransferDivMean.buoyant
#writeCellDataxyz wVarTransferDivVar.stable
#writeCellDataxyz wVarTransferDivVar.buoyant
#writeCellDataxyz wVarTransferDivMeanVar.stable
#writeCellDataxyz wVarTransferDivMeanVar.buoyant
#writeCellDataxyz wVarTransferDivZero.stable
#writeCellDataxyz wVarTransferDivZero.buoyant
#writeCellDataxyz wVarTransferSigmaDiffMean.stable
#writeCellDataxyz wVarTransferSigmaDiffMean.buoyant
#writeCellDataxyz wVarTransferSigmaDiffVar.stable
#writeCellDataxyz wVarTransferSigmaDiffVar.buoyant
#writeCellDataxyz wVarTransferSigmaDiffMeanVar.stable
#writeCellDataxyz wVarTransferSigmaDiffMeanVar.buoyant
#writeCellDataxyz wVarTransferSigmaDiffZero.stable
#writeCellDataxyz wVarTransferSigmaDiffZero.buoyant
#writeCellDataxyz wVarTransferThetaDiffMean.stable
#writeCellDataxyz wVarTransferThetaDiffMean.buoyant
#writeCellDataxyz wVarTransferThetaDiffVar.stable
#writeCellDataxyz wVarTransferThetaDiffVar.buoyant
#writeCellDataxyz wVarTransferThetaDiffMeanVar.stable
#writeCellDataxyz wVarTransferThetaDiffMeanVar.buoyant
#writeCellDataxyz wVarTransferThetaDiffZero.stable
#writeCellDataxyz wVarTransferThetaDiffZero.buoyant

writeCellDataxyz wTransferGradDiv.stable
writeCellDataxyz wTransferGradDiv.buoyant
writeCellDataxyz wTransferDrag.stable
writeCellDataxyz wTransferDrag.buoyant

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


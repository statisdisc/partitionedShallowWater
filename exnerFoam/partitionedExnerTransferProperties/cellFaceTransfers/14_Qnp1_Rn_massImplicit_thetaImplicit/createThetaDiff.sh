#!/bin/bash -e

for time in [0-9] [0-9]?? [0-9]???; do
    sumFields $time thetaDiff $time theta.stable ../00_oneFluid/$time theta.stable -scale1 -1
done


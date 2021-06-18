#!/bin/bash

# Run this file from the build directory with `source $HOME/g4PSI/build.sh`

cmake -DGeant4_DIR=$HOME/geant4_10_02_p02-install/lib64/Geant4-10.2.2/ $HOME/g4PSI
make -j8
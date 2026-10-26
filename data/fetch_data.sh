#!/bin/bash
set -e

echo "Downloading Asimov simulation data..."

wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_8.hdf5
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_16.hdf5
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_32.hdf5
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_64.hdf5
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_128.hdf5

wget -O CloudyData_UVB_HM2012.h5 "https://raw.githubusercontent.com/grackle-project/grackle_data_files/928696482fbe15d9bac4382de6134d95568f099c/input/CloudyData_UVB%3DHM2012.h5"

echo "Data download complete."

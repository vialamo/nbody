#!/bin/bash
set -e

echo "Downloading simulation data..."

# glassCube files
for size in 8 16 32 64 128; do
    FILE="glassCube_${size}.hdf5"
    echo "Attempting to download $FILE..."
    
    wget "https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/${FILE}" || \
    wget "http://www.uv.es/vicafe2/data/${FILE}"
done

# CloudyData
echo "Attempting to download CloudyData_UVB_HM2012.h5..."
wget -O CloudyData_UVB_HM2012.h5 "https://raw.githubusercontent.com/grackle-project/grackle_data_files/928696482fbe15d9bac4382de6134d95568f099c/input/CloudyData_UVB%3DHM2012.h5" || \
wget -O CloudyData_UVB_HM2012.h5 "http://www.uv.es/vicafe2/data/CloudyData_UVB_HM2012.h5"

echo "Data download complete."

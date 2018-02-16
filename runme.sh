#!/bin/bash
#
# Create all output results
#

# Useful shell settings:

# abort the script if a command fails
set -e

# abort the script if an unitialized shell variable is used
set -u

# make sure the code is up to date

pushd src
make
popd

# generate the result pictures

src/imgpro input/dogmayor.jpg output/dogmayor_0.5.jpg \
    -brightness 0.5

src/imgpro input/dogmayor.jpg output/dogmayor_1.0.jpg \
    -brightness 1.0

src/imgpro input/dogmayor.jpg output/dogmayor_1.5.jpg \
    -brightness 1.5

src/imgpro input/dogmayor.jpg output/dogmayor_sat_1.5.jpg \
        -saturation 1.5

src/imgpro input/dogmayor.jpg output/dogmayor_sat_2.0.jpg \
                -saturation 2.0

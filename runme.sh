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

# src/imgpro input/testpattern.jpg output/testpattern_harriscorners.jpg \
#                   -harris 2.0
#
src/imgpro input/test_A_alexa.jpg output/test_A_alexa_harris.jpg \
                        -harris 2.0

src/imgpro input/test_B_sony.jpg output/test_B_sony_harris.jpg \
          -harris 2.0


# src/imgpro input/test_B_sony.jpg output/test_B_sony_harris2.jpg \
#           -harris 2.0

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
# src/imgpro input/test_A_alexa.jpg output/test_A_alexa_harris.jpg \
#                         -harris 2.0
#
# src/imgpro input/test_B_sony.jpg output/test_B_sony_harris.jpg \
#           -harris 2.0
#
# src/imgpro input/testpattern.jpg output/testpattern_harris.jpg \
#           -harris 2.0

# src/imgpro input/test_E_sitting01.jpg output/test_E_sitting_match3.jpg \
#           -matchTranslation input/test_E_sitting02.jpg
#
# src/imgpro input/test_D_face01.jpg output/test_D_match3.jpg \
#           -matchTranslation input/test_D_face02.jpg

# src/imgpro input/test_C_bridge01.jpg output/test_C_harris.jpg \
#           -harris 2.0
#
# src/imgpro input/colorA.jpg output/test_C_harris2.jpg \
#           -harris 2.0


# src/imgpro input/test_C_bridge01.jpg output/test_C_bridge_match3.jpg \
#           -matchTranslation input/test_C_bridge02.jpg
#
# src/imgpro input/scc01.jpg output/scc_match.jpg \
#           -matchTranslation input/scc02.jpg


# src/imgpro input/test_B_sony.jpg output/test_B_sony_harris2.jpg \
#           -harris 2.0


# src/imgpro input/testpattern.jpg output/fake.jpg \
#           -svdTest

# src/imgpro input/colorB.jpg output/color_homography2.jpg \
#           -matchHomography input/colorA.jpg
#
# src/imgpro input/colorA.jpg output/color_homography.jpg \
#           -matchHomography input/colorB.jpg

src/imgpro input/colorA.jpg output/color_warp.jpg \
          -matchHomography input/colorB.jpg


# src/imgpro input/test_C_bridge01.jpg output/bridge_homography.jpg \
#           -matchHomography input/test_C_bridge02.jpg
# src/imgpro input/test_D_face01.jpg output/face_homography.jpg \
#           -matchHomography input/test_D_face02.jpg

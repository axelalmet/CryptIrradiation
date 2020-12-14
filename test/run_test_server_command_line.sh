#!/bin/bash

# Select which GPU to use
CUDA_VISIBLE_DEVICES=$1
# Move to the correct directory
cd ~/BuildChaste/

# Make the test
make -j$2 $3

# Run the executable with the command lines
projects/$4/test/$3 -s $5 -e $6
# # Run the test
# ctest -j$2 -R $3

# # Move the output to a non-temporary directory
# cp -r /tmp/axela/testoutput/$4 $5

# # Move LastTestLog to a non-temporary directory (this is useful for debugging)
# cp -r ~/BuildChaste/Testing/Temporary/LastTest.log $5
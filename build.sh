#!/bin/bash 


# Generate the makefiles for the project and puts them in ./build
# Then compiles and links the project with "make". The executable (if there is one) is saved in ./build/app/
# Run this with "bash build.sh"

cmake . -B build 
cd build 
make


# to use glob, cmake v3.12 or later is required (though it isn't used rn)
cmake_minimum_required(VERSION 3.12)


project(SpinChain
		VERSION 0.1) # version of my project

# check for required package Eigen3
find_package(Eigen3	REQUIRED)
		

# add subdirectories for src code and executable (app)  
# reads CMakeLists.txt file in each directory
add_subdirectory(src)
add_subdirectory(app)


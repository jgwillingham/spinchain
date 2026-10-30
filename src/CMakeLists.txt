
# Indicate the library Spin.cpp and tell it to look for includes in the ../include directory
add_library(spin Spin.cpp)
target_include_directories(spin PUBLIC ../include)

# Indicate the library SpinChain.cpp and tell it to look for includes in the ../include directory
add_library(spinchain SpinChain.cpp)
target_include_directories(spinchain PUBLIC ../include)
target_link_libraries(spinchain PUBLIC spin)




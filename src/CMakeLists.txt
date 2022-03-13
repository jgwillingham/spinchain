
# Indicate the library SpinHalf.cpp and tell it to look for includes in the ../include directory
add_library(spinhalf SpinHalf.cpp)
target_include_directories(spinhalf PUBLIC ../include)

# Indicate the library SpinHalfChain.cpp and tell it to look for includes in the ../include directory
add_library(spinhalfchain SpinHalfChain.cpp)
target_include_directories(spinhalfchain PUBLIC ../include)




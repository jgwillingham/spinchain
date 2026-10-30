
# Indicate the library sources and tell it to look for includes in the ../include directory
add_library(spinchain
    Spin.cpp
    SpinChain.cpp
    Solver.cpp
    ExactDiagonalization.cpp
    SchmidtSplit.cpp
    MPS.cpp
    MPO.cpp
    Environments.cpp
    DMRG.cpp
)
target_include_directories(spinchain PUBLIC ../include)
target_compile_features(spinchain PUBLIC cxx_std_17)

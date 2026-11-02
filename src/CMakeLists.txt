
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

# the sources include <eigen3/Eigen/...>, so the parent of Eigen's include directory is needed
target_link_libraries(spinchain PUBLIC Eigen3::Eigen)
target_include_directories(spinchain PUBLIC ${EIGEN3_INCLUDE_DIR}/..)

# the library gets linked into the shared Python extension
set_target_properties(spinchain PROPERTIES POSITION_INDEPENDENT_CODE ON)

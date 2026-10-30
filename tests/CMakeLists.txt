

### ground state tests for dimer chains

add_executable(sx_dimer_gs sx_dimer_gs.cpp)
target_include_directories(sx_dimer_gs PUBLIC ../include)
target_link_libraries(sx_dimer_gs PRIVATE spinchain)
add_test(NAME sx_dimer_gs COMMAND sx_dimer_gs)

add_executable(sy_dimer_gs sy_dimer_gs.cpp)
target_include_directories(sy_dimer_gs PUBLIC ../include)
target_link_libraries(sy_dimer_gs PRIVATE spinchain)
add_test(NAME sy_dimer_gs COMMAND sy_dimer_gs)

add_executable(sz_dimer_gs sz_dimer_gs.cpp)
target_include_directories(sz_dimer_gs PUBLIC ../include)
target_link_libraries(sz_dimer_gs PRIVATE spinchain)
add_test(NAME sz_dimer_gs COMMAND sz_dimer_gs)


### hamiltonian verification for dimer chains

add_executable(sx_dimer_hmltn sx_dimer_hmltn.cpp)
target_include_directories(sx_dimer_hmltn PUBLIC ../include)
target_link_libraries(sx_dimer_hmltn PRIVATE spinchain)
add_test(NAME sx_dimer_hmltn COMMAND sx_dimer_hmltn)

add_executable(sy_dimer_hmltn sy_dimer_hmltn.cpp)
target_include_directories(sy_dimer_hmltn PUBLIC ../include)
target_link_libraries(sy_dimer_hmltn PRIVATE spinchain)
add_test(NAME sy_dimer_hmltn COMMAND sy_dimer_hmltn)

add_executable(sz_dimer_hmltn sz_dimer_hmltn.cpp)
target_include_directories(sz_dimer_hmltn PUBLIC ../include)
target_link_libraries(sz_dimer_hmltn PRIVATE spinchain)
add_test(NAME sz_dimer_hmltn COMMAND sz_dimer_hmltn)


### hamiltonian verification for N=4 chains

add_executable(sx_4mer_hmltn sx_4mer_hmltn.cpp)
target_include_directories(sx_4mer_hmltn PUBLIC ../include)
target_link_libraries(sx_4mer_hmltn PRIVATE spinchain)
add_test(NAME sx_4mer_hmltn COMMAND sx_4mer_hmltn)

add_executable(sy_4mer_hmltn sy_4mer_hmltn.cpp)
target_include_directories(sy_4mer_hmltn PUBLIC ../include)
target_link_libraries(sy_4mer_hmltn PRIVATE spinchain)
add_test(NAME sy_4mer_hmltn COMMAND sy_4mer_hmltn)

add_executable(sz_4mer_hmltn sz_4mer_hmltn.cpp)
target_include_directories(sz_4mer_hmltn PUBLIC ../include)
target_link_libraries(sz_4mer_hmltn PRIVATE spinchain)
add_test(NAME sz_4mer_hmltn COMMAND sz_4mer_hmltn)



### spin operator algebra and spin-1 chains

add_executable(spin_operators spin_operators.cpp)
target_include_directories(spin_operators PUBLIC ../include)
target_link_libraries(spin_operators PRIVATE spinchain)
add_test(NAME spin_operators COMMAND spin_operators)

add_executable(s1_heisenberg_dimer s1_heisenberg_dimer.cpp)
target_include_directories(s1_heisenberg_dimer PUBLIC ../include)
target_link_libraries(s1_heisenberg_dimer PRIVATE spinchain)
add_test(NAME s1_heisenberg_dimer COMMAND s1_heisenberg_dimer)




### external magnetic field

add_executable(zeeman_field zeeman_field.cpp)
target_include_directories(zeeman_field PUBLIC ../include)
target_link_libraries(zeeman_field PRIVATE spinchain)
add_test(NAME zeeman_field COMMAND zeeman_field)


### numerical kernels used by DMRG

add_executable(lanczos lanczos.cpp)
target_include_directories(lanczos PUBLIC ../include)
target_link_libraries(lanczos PRIVATE spinchain)
add_test(NAME lanczos COMMAND lanczos)

add_executable(schmidt_split schmidt_split.cpp)
target_include_directories(schmidt_split PUBLIC ../include)
target_link_libraries(schmidt_split PRIVATE spinchain)
add_test(NAME schmidt_split COMMAND schmidt_split)


### matrix product states and operators

add_executable(mpo_hamiltonian mpo_hamiltonian.cpp)
target_include_directories(mpo_hamiltonian PUBLIC ../include)
target_link_libraries(mpo_hamiltonian PRIVATE spinchain)
add_test(NAME mpo_hamiltonian COMMAND mpo_hamiltonian)

add_executable(mps_state mps_state.cpp)
target_include_directories(mps_state PUBLIC ../include)
target_link_libraries(mps_state PRIVATE spinchain)
add_test(NAME mps_state COMMAND mps_state)

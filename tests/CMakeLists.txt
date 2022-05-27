

### ground state tests for dimer chains

add_executable(sx_dimer_gs sx_dimer_gs.cpp)
target_include_directories(sx_dimer_gs PUBLIC ../include)
target_link_libraries(sx_dimer_gs PRIVATE spinhalfchain)
add_test(NAME sx_dimer_gs COMMAND sx_dimer_gs)

add_executable(sy_dimer_gs sy_dimer_gs.cpp)
target_include_directories(sy_dimer_gs PUBLIC ../include)
target_link_libraries(sy_dimer_gs PRIVATE spinhalfchain)
add_test(NAME sy_dimer_gs COMMAND sy_dimer_gs)

add_executable(sz_dimer_gs sz_dimer_gs.cpp)
target_include_directories(sz_dimer_gs PUBLIC ../include)
target_link_libraries(sz_dimer_gs PRIVATE spinhalfchain)
add_test(NAME sz_dimer_gs COMMAND sz_dimer_gs)


### hamiltonian verification for dimer chains

add_executable(sx_dimer_hmltn sx_dimer_hmltn.cpp)
target_include_directories(sx_dimer_hmltn PUBLIC ../include)
target_link_libraries(sx_dimer_hmltn PRIVATE spinhalfchain)
add_test(NAME sx_dimer_hmltn COMMAND sx_dimer_hmltn)

add_executable(sy_dimer_hmltn sy_dimer_hmltn.cpp)
target_include_directories(sy_dimer_hmltn PUBLIC ../include)
target_link_libraries(sy_dimer_hmltn PRIVATE spinhalfchain)
add_test(NAME sy_dimer_hmltn COMMAND sy_dimer_hmltn)

add_executable(sz_dimer_hmltn sz_dimer_hmltn.cpp)
target_include_directories(sz_dimer_hmltn PUBLIC ../include)
target_link_libraries(sz_dimer_hmltn PRIVATE spinhalfchain)
add_test(NAME sz_dimer_hmltn COMMAND sz_dimer_hmltn)


### hamiltonian verification for N=4 chains

add_executable(sx_4mer_hmltn sx_4mer_hmltn.cpp)
target_include_directories(sx_4mer_hmltn PUBLIC ../include)
target_link_libraries(sx_4mer_hmltn PRIVATE spinhalfchain)
add_test(NAME sx_4mer_hmltn COMMAND sx_4mer_hmltn)

add_executable(sy_4mer_hmltn sy_4mer_hmltn.cpp)
target_include_directories(sy_4mer_hmltn PUBLIC ../include)
target_link_libraries(sy_4mer_hmltn PRIVATE spinhalfchain)
add_test(NAME sy_4mer_hmltn COMMAND sy_4mer_hmltn)

add_executable(sz_4mer_hmltn sz_4mer_hmltn.cpp)
target_include_directories(sz_4mer_hmltn PUBLIC ../include)
target_link_libraries(sz_4mer_hmltn PRIVATE spinhalfchain)
add_test(NAME sz_4mer_hmltn COMMAND sz_4mer_hmltn)




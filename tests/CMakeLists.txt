
add_executable(x_only_gs shc_len2_x_only.cpp)
target_include_directories(x_only_gs PUBLIC ../include)
target_link_libraries(x_only_gs PRIVATE spinhalfchain)
add_test(NAME x_only_gs_test COMMAND x_only_gs)

add_executable(y_only_gs shc_len2_y_only.cpp)
target_include_directories(y_only_gs PUBLIC ../include)
target_link_libraries(y_only_gs PRIVATE spinhalfchain)
add_test(NAME y_only_gs_test COMMAND y_only_gs)

add_executable(z_only_gs shc_len2_z_only.cpp)
target_include_directories(z_only_gs PUBLIC ../include)
target_link_libraries(z_only_gs PRIVATE spinhalfchain)
add_test(NAME z_only_gs_test COMMAND z_only_gs)



add_executable(x_only_ham shc_len2_xhamiltonian.cpp)
target_include_directories(x_only_ham PUBLIC ../include)
target_link_libraries(x_only_ham PRIVATE spinhalfchain)
add_test(NAME x_only_hamiltonian_test COMMAND x_only_ham)

add_executable(y_only_ham shc_len2_yhamiltonian.cpp)
target_include_directories(y_only_ham PUBLIC ../include)
target_link_libraries(y_only_ham PRIVATE spinhalfchain)
add_test(NAME y_only_hamiltonian_test COMMAND y_only_ham)

add_executable(z_only_ham shc_len2_zhamiltonian.cpp)
target_include_directories(z_only_ham PUBLIC ../include)
target_link_libraries(z_only_ham PRIVATE spinhalfchain)
add_test(NAME z_only_hamiltonian_test COMMAND z_only_ham)

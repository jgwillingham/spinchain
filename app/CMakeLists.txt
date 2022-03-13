

# make executable from main.cpp
add_executable(main main.cpp)

# link with spinhalfchain library (compiled SpinHalfChain.cpp)
target_link_libraries(main PRIVATE spinhalf)
target_link_libraries(main PRIVATE spinhalfchain)

target_include_directories(main PUBLIC ../include)


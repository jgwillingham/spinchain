


# make executable from main.cpp
add_executable(main main.cpp)

# link with spinchain library 
target_link_libraries(main PRIVATE spinchain)

target_include_directories(main PUBLIC ../include)


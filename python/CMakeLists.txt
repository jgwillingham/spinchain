
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core bindings.cpp)
target_link_libraries(_core PRIVATE spinchain)

install(TARGETS _core LIBRARY DESTINATION spinchain)

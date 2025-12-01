cmake_minimum_required(VERSION 3.4...3.18)
project(partial_dag_model)

add_subdirectory(third_party/pybind11)
pybind11_add_module(partial_dag_model csrc/main.cpp csrc/lasp.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(partial_dag_model
                           PRIVATE VERSION_INFO=${PDAG_VERSION_INFO})
                         
set_target_properties(partial_dag_model PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
    )                         

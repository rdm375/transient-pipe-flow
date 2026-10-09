# Native C++20 isothermal transient pipeline simulator.
#
# This module does not require or link against Fortran.
# Include it from the repository's top-level CMakeLists.txt.

option(PIPE_SIM_BUILD_CPP
       "Build the native C++20 transient simulator"
       ON)

if(PIPE_SIM_BUILD_CPP)

    add_library(pipe_sim_cpp_core STATIC
        src/cpp/physics.cpp
        src/cpp/residual.cpp
        src/cpp/jacobian.cpp
        src/cpp/banded_solver.cpp
        src/cpp/newton.cpp
        src/cpp/transient_step.cpp
        src/cpp/integrate.cpp
    )

    target_compile_features(pipe_sim_cpp_core
        PUBLIC cxx_std_20
    )

    target_include_directories(pipe_sim_cpp_core
        PUBLIC
            "${PROJECT_SOURCE_DIR}/include"
    )

    set_target_properties(pipe_sim_cpp_core PROPERTIES
        CXX_EXTENSIONS OFF
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(pipe_sim_cpp_core PRIVATE
            -Wall
            -Wextra
            -Wpedantic
        )
    endif()

    add_executable(pipe_sim_cpp
        apps/pipe_sim_cpp.cpp
    )

    target_link_libraries(pipe_sim_cpp
        PRIVATE pipe_sim_cpp_core
    )

    if(BUILD_TESTING)
        add_test(
            NAME cpp_short_transient
            COMMAND pipe_sim_cpp 100 80
        )

        add_test(
            NAME cpp_long_transient
            COMMAND pipe_sim_cpp 100 8000
        )

        add_test(
            NAME cpp_zero_steps
            COMMAND pipe_sim_cpp 100 0
        )
    endif()

endif()

# Numerical regression tests against validated reference outputs.
if(PIPE_SIM_BUILD_CPP AND BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    foreach(steps IN ITEMS 80 8000)
        add_test(
            NAME cpp_numerical_reference_${steps}
            COMMAND
                "${Python3_EXECUTABLE}"
                "${PROJECT_SOURCE_DIR}/tests/regression/check_cpp_reference.py"
                "$<TARGET_FILE:pipe_sim_cpp>"
                "${steps}"
        )
    endforeach()
endif()

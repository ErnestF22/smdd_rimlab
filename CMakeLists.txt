CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
PROJECT(degenerate_subspace_detection LANGUAGES CXX CUDA C)

set(CMAKE_CXX_FLAGS='-std=c++17')
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# add_compile_definitions(EIGEN_MALLOC_ALREADY_ALIGNED)
# add_compile_definitions(EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT)
# add_compile_definitions(EIGEN_DONT_VECTORIZE)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_86,code=sm_86)

# include(CTest)
# find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
find_package(VTK)
set(VTK_INCLUDE_DIRS /usr/include/vtk-9.1)
message("VTK_FOUND: ${VTK_FOUND}")
message("VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")
message("VTK_LIBRARIES: ${VTK_LIBRARIES}")

find_package(PCL 1.12 REQUIRED COMPONENTS io common search features visualization)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

message("PCL_FOUND: ${PCL_FOUND}")
message("PCL_INCLUDE_DIRS: ${PCL_INCLUDE_DIRS}")
message("PCL_LIBRARY_DIRS: ${PCL_LIBRARY_DIRS}")
message("PCL_DEFINITIONS: ${PCL_DEFINITIONS}")

# Require we have pkgconfig installed
find_package(PkgConfig REQUIRED)

# Tell pkgconfig to look for CSM
pkg_check_modules(CSM QUIET csm)

IF(${CSM_FOUND})
    MESSAGE(" CSM_LIBRARY_DIRS: ${CSM_LIBRARY_DIRS}")
    MESSAGE(" CSM_LIBRARIES: ${CSM_LIBRARIES}")
    MESSAGE(" CSM_INCLUDE_DIRS: ${CSM_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES(${CSM_INCLUDE_DIRS}) # important !
    LINK_DIRECTORIES(${CSM_LIBRARY_DIRS})

    # Add executable (and link with libCsm) only when libCsm exists.
    add_executable(test_fisher_matrix test/test_fisher_matrix.cpp)
    target_link_libraries(test_fisher_matrix rofl_common boost_filesystem boost_system gsl ${CSM_LIBRARIES})
    add_executable(test_fisher_matrix_gsl_eigen test/test_fisher_matrix_gsl_eigen.cpp)
    target_link_libraries(test_fisher_matrix_gsl_eigen rofl_common boost_filesystem boost_system gsl ${CSM_LIBRARIES})

# important !
ELSE(${CSM_FOUND})
    MESSAGE(ERROR " CSM not found . Check that the environment \
    variable PKG_CONFIG_PATH includes the path containing the file ’ csm.pc ’.")
ENDIF(${CSM_FOUND})

set(buildLibCsm ${CSM_FOUND})

find_package(rofl REQUIRED)
message(STATUS "rofl_FOUND ${rofl_FOUND}")
message(STATUS "rofl_INCLUDE_DIRS ${rofl_INCLUDE_DIRS}")
message(STATUS "rofl_LIBRARY_DIRS ${rofl_LIBRARY_DIRS}")
message(STATUS "rofl_LIBRARIES ${rofl_LIBRARIES}")

find_package(Eigen3 3.4 REQUIRED)
message(STATUS "Eigen3_INCLUDE_DIR: " ${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(VTK REQUIRED)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})

if(${rofl_FOUND})
    include_directories(${rofl_INCLUDE_DIRS})
    link_directories(${rofl_LIBRARY_DIRS})
endif()

include_directories(include)
include_directories(src)

add_library(gmm_estim_lib
    src/gme_MortonSort.cpp
    src/gme_GaussianMixtureEstimator.cpp
    src/dsd_utils.cpp
    src/dsd_utils_3d.cpp
    src/gme_gaussian_metric.cpp
)

# apps
add_executable(fim_ise_rtc_3d_results_interpreter app/fim_ise_rtc_3d_results_interpreter.cpp)
target_link_libraries(fim_ise_rtc_3d_results_interpreter rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(fvm3d_results_interpreter app/fvm3d_results_interpreter.cpp)
target_link_libraries(fvm3d_results_interpreter rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(run_fim_2d app/run_fim_2d.cpp)
target_link_libraries(run_fim_2d rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(run_fim_3d app/run_fim_3d.cpp)
target_link_libraries(run_fim_3d rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(run_fvm_cuda_full app/run_fvm_cuda_full.cu)
target_link_libraries(run_fvm_cuda_full rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

# target_compile_features(run_fvm_cuda_full PUBLIC cxx_std_17)
set_target_properties(run_fvm_cuda_full PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(run_fvm_cuda_full PROPERTIES CUDA_ARCHITECTURES "86")

add_executable(run_fvm app/run_fvm.cpp)
target_link_libraries(run_fvm rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(run_hessian_ise_rtc app/run_hessian_ise_rtc.cpp)
target_link_libraries(run_hessian_ise_rtc rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(run_hessian_ise_rtc_3d app/run_hessian_ise_rtc_3d.cpp)
target_link_libraries(run_hessian_ise_rtc_3d rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(run_von_mises_mixture app/run_von_mises_mixture.cpp)
target_link_libraries(run_von_mises_mixture rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

# tests

add_executable(create_degenerate_scan test/create_degenerate_scan.cpp)
target_link_libraries(create_degenerate_scan rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(estimate_pair_hessian_ise_rtc test/estimate_pair_hessian_ise_rtc.cpp)
target_link_libraries(estimate_pair_hessian_ise_rtc rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(test_fvm_cuda_full test/test_fvm_cuda_full.cu)
target_link_libraries(test_fvm_cuda_full rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

# target_compile_features(test_fvm_cuda_full PUBLIC cxx_std_17)
set_target_properties(test_fvm_cuda_full PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(test_fvm_cuda_full PROPERTIES CUDA_ARCHITECTURES "86")

add_executable(test_fvm_cuda test/test_fvm_cuda.cu)
target_link_libraries(test_fvm_cuda rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

# target_compile_features(test_fvm_cuda PUBLIC cxx_std_17)
set_target_properties(test_fvm_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(test_fvm_cuda PROPERTIES CUDA_ARCHITECTURES "86")

add_executable(test_fvm test/test_fvm.cpp)
target_link_libraries(test_fvm rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(test_gaussian_estimator_3d test/test_gaussian_estimator_3d.cpp)
target_link_libraries(test_gaussian_estimator_3d gmm_estim_lib rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(test_gme test/test_gme.cpp)
target_link_libraries(test_gme gmm_estim_lib rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(test_ise_registration test/test_ise_registration.cpp)
target_link_libraries(test_ise_registration rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})

add_executable(test_von_mises_mixture test/test_von_mises_mixture.cpp)
target_link_libraries(test_von_mises_mixture rofl_common boost_filesystem boost_system gmm_estim_lib ${PCL_LIBRARIES})



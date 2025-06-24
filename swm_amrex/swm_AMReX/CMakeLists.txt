enable_language(CUDA)
if(CMAKE_VERSION VERSION_LESS 3.20)
   include(AMReX_SetupCUDA)
endif()

set(SWM_AMREX_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/swm_mini_app_utils.cpp
    swm_mini_app_utils_amrex.cpp
)

add_executable(swm_amrex ${SWM_AMREX_SOURCES})

target_sources(swm_amrex PRIVATE ${SWM_AMREX_SOURCES})

# Or, if you have a helper function, uncomment this:
#find_package(AMReX REQUIRED CONFIG)
#include("${AMReX_DIR}/../../../AMReXCMakeModules/AMReXTargetHelpers.cmake")
#include("/glade/u/home/htorres/amrex_build_GNU_MPI_CUDA/install/lib/cmake/AMReX/AMReXCMakeModules/AMReXTargetHelpers.cmake")
include("AMReXTargetHelpers")
setup_target_for_cuda_compilation(swm_amrex)

target_include_directories(swm_amrex PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex)
target_include_directories(swm_amrex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(swm_amrex PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex PUBLIC m)

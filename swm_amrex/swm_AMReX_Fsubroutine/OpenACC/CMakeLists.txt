#find_package(AMReX REQUIRED FORTRAN) 

# TODO: This is basicaly a copy of the block from the CMakeLists.txt for the C OpenACC build. Maybe combine in a function?
# Set OpenACC flags based on compiler and device
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc -foffload=disable)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -fopenacc)
    endif()
    set(OPENACC_LINK_FLAGS -fopenacc)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=multicore)
    elseif(SWM_DEVICE STREQUAL "gpu")
        #set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma)
        set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma -Mbounds -gpu=debug)
        #set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma -mcmodel=medium)
        set(OPENACC_LINK_FLAGS -acc=gpu)
    endif()
else()
    # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
    message(FATAL_ERROR "We have not added support for your compiler (${CMAKE_Fortran_COMPILER_ID}) for OpenACC build. Exiting configuration.")
    #message(WARNING "We have not added your support for your compiler (${CMAKE_Fortran_COMPILER_ID}) for OpenACC build. So we are not passing any OpenACC flags to the compiler. It is likely the pragmas will be ignored.")
    #set(OPENACC_COMPILE_FLAGS "")
    #set(OPENACC_LINK_FLAGS "")
endif()

################################################################################
## Fortran compiled seperatly 
################################################################################

set(SWM_AMREX_FSUBROUTINE_ACC_CPP_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/swm_mini_app_utils.cpp
    ../swm_mini_app_utils_fsubroutine.cpp
)

add_library(swm_fortran_kernels_obj OBJECT swm_fortran_kernelsACC.F90)
# Maybe also add this flag... -Mnomain
set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_OPTIONS ${OPENACC_COMPILE_FLAGS})
#set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_FLAGS "${OPENACC_COMPILE_FLAGS}")
set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES
    LANGUAGE Fortran
    CUDA_SOURCE_PROPERTY_FORMAT OFF
)
target_include_directories(swm_fortran_kernels_obj PUBLIC $<TARGET_PROPERTY:AMReX::amrex_2d,INTERFACE_INCLUDE_DIRECTORIES>)

if (SWM_CUDA)
  # Prevent CMake from passing NVCC arch flags to nvc++
  set(CMAKE_CUDA_ARCHITECTURES "" CACHE STRING "No CUDA arch flags for nvc++" FORCE)
  set(CUDA_ARCHITECTURES "" CACHE STRING "No CUDA arch flags for nvc++" FORCE)
  # Set CUDA host linker to nvc++ so -acc=gpu is understood
  set(CMAKE_CUDA_HOST_LINK_LAUNCHER "/glade/u/apps/casper/24.12/spack/opt/spack/nvhpc/24.11/gcc/12.4.0/bzc3/Linux_x86_64/24.11/compilers/bin/nvc++")
  set(CMAKE_CUDA_LINK_EXECUTABLE "/glade/u/apps/casper/24.12/spack/opt/spack/nvhpc/24.11/gcc/12.4.0/bzc3/Linux_x86_64/24.11/compilers/bin/nvc++ <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
endif()

add_executable(swm_amrex_fsubroutine_acc
    ${SWM_AMREX_FSUBROUTINE_ACC_CPP_SOURCES}
    $<TARGET_OBJECTS:swm_fortran_kernels_obj>
)

# Now that the target exists, set properties and link options
if (SWM_CUDA)
  # Remove NVCC arch flags and add -gpu=cc80 for nvc++
  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES CUDA_ARCHITECTURES "")
  target_link_options(swm_amrex_fsubroutine_acc PRIVATE -gpu=cc80)
  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES LINKER_LANGUAGE CUDA)
else()
  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES LINKER_LANGUAGE CXX)
endif()

target_link_options(swm_amrex_fsubroutine_acc PRIVATE ${OPENACC_LINK_FLAGS})

if (SWM_CUDA)
  # Ensure no CUDA arch flags are re-added by AMReX helpers
  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES CUDA_ARCHITECTURES "")
  include("AMReXTargetHelpers")
  setup_target_for_cuda_compilation(swm_amrex_fsubroutine_acc
    SOURCES ${SWM_AMREX_FSUBROUTINE_ACC_CPP_SOURCES}
  )
  # Clear again in case AMReX helpers re-add it
  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES CUDA_ARCHITECTURES "")
endif()

target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex)
target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC m)
#target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC nvf)

###############################################################################
# Everything together
###############################################################################

#set(CMAKE_CUDA_HOST_LINK_LAUNCHER "/glade/u/apps/casper/24.12/spack/opt/spack/nvhpc/24.11/gcc/12.4.0/bzc3/Linux_x86_64/24.11/compilers/bin/nvc++")
#
#set(SWM_AMREX_FSUBROUTINE_ACC_SOURCES
#    ${PROJECT_SOURCE_DIR}/swm_amrex/main.cpp
#    ${PROJECT_SOURCE_DIR}/swm_amrex/swm_mini_app_utils.cpp
#    ../swm_mini_app_utils_fsubroutine.cpp
#    swm_fortran_kernelsACC.F90
#)
#
#
#add_executable(swm_amrex_fsubroutine_acc
#    ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES}
#)
#
#set(CMAKE_CXX_LINK_EXECUTABLE "/glade/u/apps/casper/24.12/spack/opt/spack/nvhpc/24.11/gcc/12.4.0/bzc3/Linux_x86_64/24.11/compilers/bin/nvc++ <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
#
#set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES LANGUAGE CUDA)
#
#if(SWM_CUDA)
#  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES LINKER_LANGUAGE CUDA)
#else()
#  set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES LINKER_LANGUAGE CXX)
#endif()
##set_target_properties(swm_amrex_fsubroutine_acc PROPERTIES LINKER_LANGUAGE CXX)
#
#set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_OPTIONS "${OPENACC_COMPILE_FLAGS}")
#
#
#target_link_options(swm_amrex_fsubroutine_acc PRIVATE ${OPENACC_LINK_FLAGS})
##target_link_options(swm_amrex_fsubroutine_acc
##  PRIVATE
##    $<$<OR:$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,Fortran>,$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,CUDA>>:${OPENACC_LINK_FLAGS}>
##)
##if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
##  target_link_options(swm_amrex_fsubroutine_acc PRIVATE ${OPENACC_LINK_FLAGS})
##endif()
#
#if (SWM_CUDA)
#  include("AMReXTargetHelpers")
#  setup_target_for_cuda_compilation(swm_amrex_fsubroutine_acc
#    SOURCES ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES}
#  )
#endif()
#
#target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex)
#target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
#
#target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC AMReX::amrex_2d)
#target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC m)
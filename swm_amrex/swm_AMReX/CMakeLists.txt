set(SWM_AMREX_SOURCES
    main.cpp
    swm_mini_app_kernels.h
    swm_mini_app_utils.cpp
    swm_mini_app_utils.h
)

add_executable(swm_amrex ${SWM_AMREX_SOURCES})

target_link_libraries(swm_amrex PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex PUBLIC m)
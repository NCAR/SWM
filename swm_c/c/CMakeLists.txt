###############################################################################
# Plain C version
###############################################################################
if(SWM_DEVICE STREQUAL "cpu")

  set(SWM_C_SOURCES
      shallow_swap.c
      wtime.c
  )
  
  add_executable(swm_c ${SWM_C_SOURCES})
  
  # Link math library
  target_link_libraries(swm_c m)

endif()

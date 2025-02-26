#ifndef AMREX_HAVE_NO_CONFIG_H
#define AMREX_HAVE_NO_CONFIG_H
#ifndef AMREX_GPU_MAX_THREADS
#define AMREX_GPU_MAX_THREADS 0
#endif
#ifndef BL_SPACEDIM
#define BL_SPACEDIM 2
#endif
#ifndef AMREX_SPACEDIM
#define AMREX_SPACEDIM 2
#endif
#ifndef BL_FORT_USE_UNDERSCORE
#define BL_FORT_USE_UNDERSCORE 1
#endif
#ifndef AMREX_FORT_USE_UNDERSCORE
#define AMREX_FORT_USE_UNDERSCORE 1
#endif
#ifndef BL_Linux
#define BL_Linux 1
#endif
#ifndef AMREX_Linux
#define AMREX_Linux 1
#endif
#ifndef NDEBUG
#define NDEBUG 1
#endif
#ifdef __cplusplus
#endif
#if defined(AMREX_USE_OMP) && !defined(_OPENMP) && !defined(AMREX_USE_HIP)
#error libamrex was built with OpenMP
#endif
#endif

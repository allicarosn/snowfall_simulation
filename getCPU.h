///// header file for CPU timing function /////
///// now compatible with openMP and MPI  /////

#ifndef GET_CPU_H
#define GET_CPU_H

#include <ctime>
// ---------------------------------------------
// Return the current wall-clock time in seconds
// ---------------------------------------------
inline double getCPU()
{
  #ifdef _OPENMP
  return omp_get_wtime(); // use omp timer
  #elifdef USE_PPP
    return MPI_Wtime(); // use MPI timer
  #else
    return ( 1.0*std::clock() )/CLOCKS_PER_SEC ; // usual timer
  #endif
}

#endif


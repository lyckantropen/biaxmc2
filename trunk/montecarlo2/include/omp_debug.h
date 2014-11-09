#ifndef __OMP_DEBUG
#define __OMP_DEBUG

/// ONLY FOR DEBUGGING
/// on Mac OS X clang does not support OpenMP
/// so we provide dummy functions
/// so that the code compiles
/// but no parallelism

int omp_get_num_threads();
int omp_get_thread_num();
void omp_set_dynamic(int);
void omp_set_num_threads(int);

#endif

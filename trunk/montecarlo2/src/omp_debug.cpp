#ifdef OMP_DEBUG
#include "omp_debug.h"

int omp_get_num_threads()
{
    return 1;
}


int omp_get_thread_num()
{
    return 0;
}


void omp_set_dynamic(int)
{
    //
}


void omp_set_num_threads(int)
{
    //
}

#endif

/* 
 * File:   main.cpp
 * Author: karol
 *
 * Created on 4 listopad 2009, 17:11
 */

#include "boost.h"
#include "std.h"
#include "exceptions.h"
#include "serialsaveload.h"
#include "quantity_iterator.h"
#include "Particle.h"
#include "Lattice.h"
#include "Metropolis.h"
#include "PRE79StandardHamiltonian.h"
#include "PRE79StandardProperties.h"
#include "Settings.h"
#include "SimulationDB.h"
#include "PRE79Scanning.h"
//#include <omp.h>

int main(int argc, char** argv)
{
    std::string cfg("simulation.ini");
    if(argc>1)
        cfg=argv[1];
    Settings setup(cfg);
    if(setup.scanning.enabled){
        PRE79Scanning scanning(setup);
        scanning.SetStream(&std::cout);
        if(setup.scanning.threaded)
            scanning.RunParallel();
        else
            scanning.RunNonParallel();
    }
    else {
        PRE79Simulation simulation(setup);
        simulation.SetStream(&std::cout);
        simulation.Run();
    }
    return (EXIT_SUCCESS);
}


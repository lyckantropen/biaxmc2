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
#include "PRE79Simulation.h"
#include "Settings.h"
#include "SimulationDB.h"
#include "PRE79Scanning.h"
//#include <omp.h>

int main(int argc, char** argv)
{
    Settings setup("/home/karol/simulation.ini");
    SimulationDB db(setup);
    if(setup.scanning.enabled){
        PRE79Scanning scanning(setup,db);
        scanning.SetStream(&std::cout);
        if(setup.scanning.threaded)
            scanning.RunParallel();
        else
            scanning.RunNonParallel();
    }
    else {
        PRE79Simulation simulation(setup,db);
        simulation.SetStream(&std::cout);
        simulation.Run();
    }
    return (EXIT_SUCCESS);
}


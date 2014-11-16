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
#include "SimulationDB.h"
#include "PRE79Scanning.h"
#include "ParallelTempering.h"
#include "cstdlib"

#ifndef OMP_DEBUG
#include <omp.h>
#else
#include "omp_debug.h"
#endif

int main(int argc, char** argv)
{
    //testowy komentarz
    setenv("OMP_SCHEDULE", "static", 0);
    std::string cfg("simulation.ini");
    std::string mode("submit");

    if(argc > 1)
        cfg = argv[1];

    if(cfg == std::string("help"))
    {
        std::cout << "Syntax: " << std::endl
                  << "To run the simulation with options set in 'simulation.ini' in place: \n\t$ " << argv[0] << " simulation.ini run " << std::endl
                  << "To execute the simulation in a PBS queue (Shiva cluster supported): \n\t$ " << argv[0] << " simulation.ini" << std::endl;

        std::cout << std::endl << "Configuration file option details: " << std::endl;
        std::cout << Settings().ListOptions();
        return(EXIT_SUCCESS);
    }

    Settings setup(cfg);
    setup.SetStream(&std::cout);

    if(argc > 2)
        mode = argv[2];

    if(mode == "run")
    {

        if(setup.scanning.enabled)
        {
            if(setup.scanning.parallel_tempering)
            {
                ParallelTempering par(setup);
                par.Run2();
            }
            else
            {
                PRE79Scanning scanning(setup);
                scanning.SetStream(&std::cout);
                if(setup.scanning.threaded)
                    if(setup.scanning.threaded_production)
                        scanning.RunParallelParallel();
                    else
                        scanning.RunParallel();
                else
                    scanning.RunNonParallel();
            }
        }
        else
        {
            PRE79Simulation simulation(setup);
            simulation.SetStream(&std::cout);
            simulation.Run();
        }
    }
    else
    {
        std::stringstream pt;
        std::string exec = argv[0];


        pt      << "#!/bin/bash\n"
                << "#PBS -N " << setup.project.name << std::endl
                << "#PBS -l cput=" << setup.pbs.cput << std::endl
                << "#PBS -q " << setup.pbs.queue << std::endl
                << "#PBS -m ae\n"
                << std::endl
                << "cd " << fs::current_path().string() << std::endl
                << exec << " " << cfg << " run >> " << setup.project.name << ".txt" << std::endl
                << std::endl;

        std::string script = pt.str();

        std::stringstream s;
        s << "pbs-" << setup.project.name << std::rand() ;
        fs::path        script_temp(s.str());
        fs::path        tmp(".tmp");
        if(!fs::exists(tmp))
            fs::create_directory(tmp);

        std::ofstream   t((tmp / script_temp).string().c_str());
        t << script ;
        t.close();

        if(mode != "generate")
        {
            std::string     sub_command("qsub ");
            sub_command += (tmp / script_temp).string();
            std::system(sub_command.c_str());
            fs::remove(tmp / script_temp);
        }

        std::cout << script;


    }
    return (EXIT_SUCCESS);
}


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
//#include <omp.h>

int main(int argc, char** argv)
{
    setenv("OMP_SCHEDULE","static",0);
    std::string cfg("simulation.ini");
    std::string mode("submit");

    if(argc>1)
        cfg=argv[1];
    Settings setup(cfg);
    setup.SetStream(&std::cout);

    //rng_wrap<rangen_mwc>::setup(setup.openmp.number_of_threads);
    
    if(argc>2)
        mode=argv[2];

    if(mode=="run"){
            
        if(setup.scanning.enabled){
            if(setup.scanning.parallel_tempering){
                ParallelTempering par(setup);
                par.Run2();
            }
            else {
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
        else {
            PRE79Simulation simulation(setup);
            simulation.SetStream(&std::cout);
            simulation.Run();
        }
    }
    else {
        std::stringstream pt;
        std::string cput;
        std::string exec=argv[0];
        if(setup.pbs.queue=="normal")
            cput="12:00:00";
        if(setup.pbs.queue=="long")
            cput="72:00:00";
        if(setup.pbs.queue=="vlong")
            cput="10000:00:00";
        if(setup.pbs.queue=="mp4")
            cput="340:00:00";
        if(setup.pbs.queue=="mp8")
            cput="2000:00:00";
        if(setup.pbs.queue=="mp16")
            cput="900:00:00";
        if(setup.pbs.queue=="mp24")
            cput="1200:00:00";

        
        pt      << "#!/bin/bash\n"
                << "#PBS -N "<< setup.project.name << std::endl
                << "#PBS -l cput=" << cput <<std::endl
                << "#PBS -q "<< setup.pbs.queue << std::endl
                << "#PBS -m ae\n"
                << std::endl
                << "cd "<<fs::current_path().string() << std::endl
                << exec << " " << cfg << " run >> "<< setup.project.name << ".txt" << std::endl
                << std::endl;

        std::string script = pt.str();

        std::stringstream s;
        s << "pbs-"<< setup.project.name << std::rand() ;
        fs::path        script_temp(s.str());
        fs::path        tmp(".tmp");
        if(!fs::exists(tmp))
                fs::create_directory(tmp);

        std::ofstream   t((tmp/script_temp).string().c_str());
        t << script ;
        t.close();

        if(mode!="generate"){
            std::string     sub_command("qsub ");
            sub_command+=(tmp/script_temp).string();
            std::system(sub_command.c_str());
            fs::remove(tmp/script_temp);
        }

        std::cout << script;


    }
    return (EXIT_SUCCESS);
}


/* 
 * File:   PRE79Scanning.h
 * Author: karol
 *
 * Created on 7 grudzień 2009, 22:09
 */

#ifndef _PRE79SCANNING_H
#define	_PRE79SCANNING_H

#include "PRE79Simulation.h"
#include "SimulationDB.h"
#include "Settings.h"
#include "ILoggable.h"
#include <omp.h>

class PRE79Scanning:public ILoggable {
    SimulationDB    &    db;                ///<baza danych
    const Settings &        settings;       ///<globalne ustawienia
    //PRE79Simulation         * simulation;   ///<bieżąca symulacja
    
    const double start;
    const double end;
    const double delta;
    const int nscans;
    const std::string variable;
public:
    PRE79Scanning(const Settings & set,SimulationDB & dbase):
    db(dbase),
    settings(set),
    start(set.scanning.start),
    end(set.scanning.end),
    delta(set.scanning.delta),
    nscans((end-start)/delta+2),
    variable(set.scanning.variable)
    {
        //simulation = new PRE79Simulation(settings,db);
    }               
    void RunNonParallel(){
        Lattice state;
        Settings current_settings = settings ;
	Log() << "Calculating expected time of execution\n";
        pt::time_duration timeof1000cycles = PRE79Simulation(settings,db).DurationOf1000Cycles();
        int nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles + settings.simulation.thermalization_cycles)/1000;
        Log() << "Expected time of execution: " << pt::to_simple_string(timeof1000cycles*nkcycles) << std::endl;
        
        Log() << "Non-parallel version\n";
        Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl; 
        
        // wersja bez paralelizacji, bo mamy ładowanie poprzedniego stanu
        //TODO: pozbyć się "simulation" i tworzyć symulację w każdej iteracji
        pt::ptime start_t = pt::second_clock::local_time();
        for(int i=0;i<nscans;i++){
            double value = start + double(i)*delta;
            Log() << "Value: " << value << std::endl;

            if(variable=="hamiltonian.tau")
                current_settings.hamiltonian.tau=value;
            if(variable=="hamiltonian.lambda")
                current_settings.hamiltonian.lambda=value;
            if(variable=="hamiltonian.temperature")
                current_settings.hamiltonian.temperature=value;


            if(i==0){
                PRE79Simulation simulation(current_settings,db);
                simulation.SetStream(&Log());
                state = *simulation.Run();
            }
            else{
                PRE79Simulation simulation(current_settings,db,state);
                simulation.SetStream(&Log());
                state = *simulation.Run();
            }

            /*
            simulation->SetStream(&Log());
            state = *simulation->Run();
            //Log() << simulation->GetLog();
            delete simulation;
            simulation = new PRE79Simulation(current_settings,db,state);
            
            double value = start + double(i+1)*delta;
            */
            //międzyczasy
            pt::time_duration t1scan = pt::second_clock::local_time()-start_t;
            std::cout << "Total " << pt::to_simple_string(t1scan*(nscans-i-1)) << " remaining\n";
            start_t = pt::second_clock::local_time();
        }
        Log() << "Scanning finished\n";
    }
    void RunParallel(){
        Log() << "Parallel OpenMP version\n";
        Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl;
        Log() << "Trying to set the number of concurrent simulations to " << settings.scanning.number_of_threads << std::endl;

	omp_set_dynamic(1);
        omp_set_num_threads(settings.scanning.number_of_threads);
        #pragma omp parallel for
        for(int i=0;i<nscans;i++){

            double value = start + double(i)*delta;

            Settings current_settings(settings);

            if(variable=="hamiltonian.tau")
                current_settings.hamiltonian.tau=value;
            if(variable=="hamiltonian.lambda")
                current_settings.hamiltonian.lambda=value;
            if(variable=="hamiltonian.temperature")
                current_settings.hamiltonian.temperature=value;

            Log() << "Thread: "<< omp_get_thread_num() << ", Value: " << value << std::endl;
            PRE79Simulation thread_sim(current_settings,db);
            thread_sim.SetStream(&Log());
            thread_sim.Run();
            
        }
        Log() << "Scanning finished\n";
    }
    ~PRE79Scanning(){
        //delete simulation;
    }
};


#endif	/* _PRE79SCANNING_H */


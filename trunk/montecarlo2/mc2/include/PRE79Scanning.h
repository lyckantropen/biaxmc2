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
#include "SimulationDBFind.h"
#include <omp.h>

class PRE79Scanning:public ILoggable {
    const Settings &        settings;       ///<globalne ustawienia
    //PRE79Simulation         * simulation;   ///<bieżąca symulacja
    
    const double start;
    const double end;
    const double delta;
    const int nscans;
    const std::string variable;
public:
    PRE79Scanning(const Settings & set):
    settings(set),
    start(set.scanning.start),
    end(set.scanning.end),
    delta(set.scanning.delta),
    nscans((end-start)/delta+2),
    variable(set.scanning.variable)
    {
    }
    /**
     * Uruchamianie symulacji dla różnych wartości parametrów. Tutaj może wejść
     * paralelizacja na poziomie produkcji!
     */
    void RunNonParallel(){
        Lattice state;
        Settings current_settings = settings ;

        //--- obliczanie czasu działania
	Log() << "Calculating expected time of execution\n";
        pt::time_duration timeof1000cycles = PRE79Simulation(settings).DurationOf1000Cycles();
        int nkcycles=0;
        if(settings.scanning.threaded_production)
            nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles/settings.openmp.number_of_threads + settings.simulation.thermalization_cycles)/1000;
        else
            nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles + settings.simulation.thermalization_cycles)/1000;

        Log() << "Expected time of execution: " << pt::to_simple_string(timeof1000cycles*nkcycles) << std::endl;
        //---


        Log() << "Non-parallel version\n";
        Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl; 
        
        // wersja bez paralelizacji, bo mamy ładowanie poprzedniego stanu
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

            //--- znajdowanie zapisanego stanu - jeżeli znajdziemy, możemy pominąć ten krok
            if(settings.scanning.continue_if_results_exist){
                Log() << "Searching database for previous state\n";
                bool found=false;
                Lattice found_state = FindFinalState(current_settings,found);
                if(found){
                    state = found_state;
                    Log() << "Previous state recovered, skipping\n";
                    // pomijamy ten krok symulacji, stan o dokładnie takich samych parametrach został już policzony
                    continue;
                }
                else{
                    Log() << "No previous state recovered, going on\n";
                }
            }
            //---


            if(i==0 || settings.scanning.reuse_thermalized==false){
                PRE79Simulation simulation(current_settings);
                Log() << simulation.GetLog();
                simulation.SetStream(&Log());
                if(settings.scanning.threaded_production)
                    state = simulation.RunParallel();
                else
                    state = simulation.Run();
            }
            else{
                PRE79Simulation simulation(current_settings,state);
                Log() << simulation.GetLog();
                simulation.SetStream(&Log());
                if(settings.scanning.threaded_production)
                    state = simulation.RunParallel();
                else
                    state = simulation.Run();
            }
            //międzyczasy
            pt::time_duration t1scan = pt::second_clock::local_time()-start_t;
            std::cout << "Total " << pt::to_simple_string(t1scan*(nscans-i-1)) << " remaining\n";
            start_t = pt::second_clock::local_time();
        }
        Log() << "Scanning finished\n";
    }
    /**
     * Naiwna paralelizacja, równoważna uruchomieniu wielu symulacji naraz.
     */
    void RunParallel(){
        Log() << "Parallel OpenMP version\n";
        Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl;
        Log() << "Trying to set the number of concurrent simulations to " << settings.openmp.number_of_threads << std::endl;

        if(settings.openmp.dynamic)
            omp_set_dynamic(1);
        else
            omp_set_dynamic(0);
        omp_set_num_threads(settings.openmp.number_of_threads);
        #pragma omp parallel for
        for(int i=0;i<nscans;i++){

            double value = start + double(i)*delta;
            #pragma omp critical
            Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ", Value: " << value << std::endl;

            Settings current_settings(settings);

            if(variable=="hamiltonian.tau")
                current_settings.hamiltonian.tau=value;
            if(variable=="hamiltonian.lambda")
                current_settings.hamiltonian.lambda=value;
            if(variable=="hamiltonian.temperature")
                current_settings.hamiltonian.temperature=value;

            //--- znajdowanie zapisanego stanu - jeżeli znajdziemy, możemy pominąć ten krok
            if(settings.scanning.continue_if_results_exist){
                #pragma omp critical
                Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": Searching database for results of previous simulation\n";
                bool found=false;
                FindFinalProperties(current_settings,found);
                if(found){
                    #pragma omp critical
                    Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": Previous results recovered, skipping\n";
                    // pomijamy ten krok symulacji, stan o dokładnie takich samych parametrach został już policzony
                    continue;
                }
                else{
                    #pragma omp critical
                    Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": No previous results recovered, going on\n";
                }
            }
            //---

            PRE79Simulation thread_sim(current_settings);
            #pragma omp critical
            Log() << thread_sim.GetLog();
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


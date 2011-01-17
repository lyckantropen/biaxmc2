/* 
 * File:   PRE79Scanning.h
 * Author: karol
 *
 * Created on 7 grudzień 2009, 22:09
 */

#ifndef _PRE79SCANNING_H
#define	_PRE79SCANNING_H

//#include "PRE79Simulation.h"
#include "Settings.h"
#include "ILoggable.h"
#include "RuntimePropertiesServer.h"
#include <omp.h>

//void call_run(Lattice & state, PRE79Simulation & sim) ;

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
        Settings current_settings(settings) ;

        //--- obliczanie czasu działania
	if(settings.simulation.calculate_time){
		Log() << "Calculating expected time of execution\n";
		pt::time_duration timeof1000cycles = PRE79Simulation(settings).DurationOf1000Cycles();
		int nkcycles=0;
		if(settings.scanning.threaded_production)
		    nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles/settings.openmp.number_of_threads + settings.simulation.thermalization_cycles)/1000;
		else
		    nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles + settings.simulation.thermalization_cycles)/1000;

		Log() << "Expected time of execution: " << pt::to_simple_string(timeof1000cycles*nkcycles) << std::endl;
	}
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
	    if(variable=="hamiltonian.h")
		current_settings.hamiltonian.h=value;

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
                //serwer dostępu do symulacji w trakcie jej trwania, uruchamiany w osobnym wątku
                if(settings.output.start_service){
                    RuntimePropertiesServer rt(current_settings,simulation);
                    boost::thread thread1(rt);
                }
                //--
                Log() << simulation.GetLog();
                simulation.SetStream(&Log());
                if(settings.scanning.threaded_production)
                    state = simulation.RunParallel();
                else
                    state = simulation.Run();
                    
            }
            else{
                PRE79Simulation simulation(current_settings,state);
                //serwer dostępu do symulacji w trakcie jej trwania, uruchamiany w osobnym wątku
                if(settings.output.start_service){
                    RuntimePropertiesServer rt(current_settings,simulation);
                    boost::thread thread1(rt);
                }
                //--
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
        #pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
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
            if(current_settings.scanning.continue_if_results_exist){
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
    /*
     * Dwustopniowa paralelizacja: najpierw termalizacje, potem produkcje
     */
    void RunParallelParallel() {
        Log() << "Parallel thermalization with parallel productions\n";
        Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl;
        Log() << "Trying to set the number of concurrent simulations to " << settings.openmp.number_of_threads << std::endl;

        if(settings.openmp.dynamic)
            omp_set_dynamic(1);
        else
            omp_set_dynamic(0);
        omp_set_num_threads(settings.openmp.number_of_threads);


        //-- niesparalelizowana pętla, jedziemy kawałkami rozmiaru liczby dostępnych wątków
        int incr=settings.openmp.number_of_threads;
        for(int i=0; i<nscans; i+=incr){

            int chunk = i/incr+1;
            int chunks = nscans/incr+1;

            int current_nscans = incr;
            if(chunk*incr > nscans)
                current_nscans = nscans-(chunk-1)*incr;
            //-- stany do stermalizowania
            std::vector<Lattice>    states(current_nscans);


            Log() << "Chunk " << chunk << "/" << chunks << ": starting\n";
            Log() << "Chunk " << chunk << "/" << chunks << ": that is " << current_nscans << " thermalizations/productions of " << nscans << " total\n";
            Log() << "Chunk " << chunk << "/" << chunks << ": thermalization\n";

            //-- termalizacja
            #pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
            for(int t=0;t<current_nscans;t++){
                double value = start + double((chunk-1)*incr+t)*delta;
                #pragma omp critical
                Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ", Value: " << value << std::endl;

                Settings current_settings(settings);

                if(variable=="hamiltonian.tau")
                    current_settings.hamiltonian.tau=value;
                if(variable=="hamiltonian.lambda")
                    current_settings.hamiltonian.lambda=value;
                if(variable=="hamiltonian.temperature")
                    current_settings.hamiltonian.temperature=value;

                //current_settings.simulation.production_cycles=1;
                //current_settings.openmp.number_of_threads=1;
                //current_settings.scanning.enabled=false;
                //current_settings.scanning.threaded=false;
                //current_settings.scanning.threaded_production=false;

                PRE79Simulation termo(current_settings);
                #pragma omp critical
                Log() << termo.GetLog();
                termo.SetStream(&Log());
                states[t]=termo.Thermalize();
            }

            Log() << "Chunk " << chunk << "/" << chunks << ": finished thermalization\n";
            Log() << "Chunk " << chunk << "/" << chunks << ": starting production\n";

            //-- produkcja
            //#//pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
            for(int t=0;t<current_nscans;t++){
                double value = start + double((chunk-1)*incr+t)*delta;
                #pragma omp critical
                Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ", Value: " << value << std::endl;

                Settings current_settings(settings);

                if(variable=="hamiltonian.tau")
                    current_settings.hamiltonian.tau=value;
                if(variable=="hamiltonian.lambda")
                    current_settings.hamiltonian.lambda=value;
                if(variable=="hamiltonian.temperature")
                    current_settings.hamiltonian.temperature=value;

                current_settings.simulation.thermalization_cycles=0;
                current_settings.simulation.supplementary_thermalization_cycles=0;


                //--- znajdowanie zapisanego stanu - jeżeli znajdziemy, możemy pominąć ten krok
                if(current_settings.scanning.continue_if_results_exist){
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


                PRE79Simulation production(current_settings,states[t]);
                #pragma omp critical
                Log() << production.GetLog();
                production.SetStream(&Log());
                production.RunParallel();
            }
            Log() << "Chunk " << chunk << "/" << chunks << ": finished\n";
        }
        Log() << "Scanning finished\n";
    }

    ~PRE79Scanning(){
        //delete simulation;
    }
};


#endif	/* _PRE79SCANNING_H */


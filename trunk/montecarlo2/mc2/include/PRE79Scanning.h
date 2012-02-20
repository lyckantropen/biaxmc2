/** 
 * @file PRE79Scanning.h
 * 
 *
 */

#ifndef _PRE79SCANNING_H
#define	_PRE79SCANNING_H

#include "Settings.h"
#include "RuntimePropertiesServer.h"
#include <omp.h>

/**
 * Manages different variants of scanning over a set of parameters. This is an
 * ugly class. Definitely some things could have been done more elegantly.
 */
class PRE79Scanning:public ILoggable {
    Settings &        settings;         ///<a local reference to the global settings object
    
    const double start;                 ///<start value
    const double end;                   ///<end value
    const double delta;                 ///<step
    const int nscans;                   ///<total number of parameter values
    const std::string variable;         ///<name of the parameter to scan over
public:
    /**
     * Constructor
     * @param set Global settings object
     */
    PRE79Scanning(Settings & set):
    settings(set),
    start(set.scanning.start),
    end(set.scanning.end),
    delta(set.scanning.delta),
    nscans((end-start)/delta+2),
    variable(set.scanning.variable)
    {
    }

    /**
     * The most trivial scanning possible. Single-threaded scanning procedure, 
     * where simulations are run one by one. A simulation is started, the state
     * is thermalized, the production is carried out and the simulation for the
     * next value of the parameter isn't started until the previous one has ended.
     * 
     * Stil, there is possibility for some parallelism, namely parallel production.
     * The only change in the picture is that the originally single-threaded
     * production is split into threads. The final state with the least energy
     * is then passed as initial state to the next simulation for the next
     * parameter value.
     */
    void RunNonParallel(){
        ///the state that will be passed on between parameter values
        Lattice state;
        ///we keep a copy of "localized" settings, where the actual present values are
        Settings current_settings(settings) ;

        ///if needed, we can estimate the time of the execution
	if(settings.simulation.calculate_time){
		Log() << "Calculating expected time of execution\n";
                ///the estimation is based on the time of execution of 1000 cycles
                ///so in some situations it can be inconvenient
		pt::time_duration timeof1000cycles = PRE79Simulation(settings).DurationOf1000Cycles();
		int nkcycles=0;
		if(settings.scanning.threaded_production)
		    nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles/settings.openmp.number_of_threads + settings.simulation.thermalization_cycles)/1000;
		else
		    nkcycles = ((nscans-1)*settings.simulation.supplementary_thermalization_cycles + nscans*settings.simulation.production_cycles + settings.simulation.thermalization_cycles)/1000;

		Log() << "Expected time of execution: " << pt::to_simple_string(timeof1000cycles*nkcycles) << std::endl;
	}



        Log() << "Non-parallel version\n";
        ///here we look if we are iterating over predefined values or over an interval with a specified step
        if(settings.scanning.separate_values)
            Log() << "Scanning values: " << settings.scanning.values << std::endl;
        else
            Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl; 
        
        
      
        ///initialize the clock for estimating the remaining time
        pt::ptime start_t = pt::second_clock::local_time();
        
        ///loop over all parameter values
        for(int i=0;i<nscans;i++){
            
            ///this is the actual present value of the parameter
            double value = 0.0;
            ///we set it to the appropriate value depending on the type of the set we are scanning over
            if(settings.scanning.separate_values)
                value = settings.scanning.values[i];
            else
                value = start + double(i)*delta;
            Log() << "Value: " << value << std::endl;

            ///several parameter names are supported
            ///we alter the "current settings" with accordance to the specified scanning variable name
            if(variable=="hamiltonian.tau")
                current_settings.hamiltonian.tau=value;
            if(variable=="hamiltonian.lambda")
                current_settings.hamiltonian.lambda=value;
            if(variable=="hamiltonian.temperature")
                current_settings.hamiltonian.temperature=value;
	    if(variable=="hamiltonian.h")
		current_settings.hamiltonian.h=value;

            ///here we search the database for an existing result for the present value
            if(settings.scanning.continue_if_results_exist){
                Log() << "Searching database for previous state\n";
                bool found=false;
                ///the actual search is implemented in the following function
                Lattice found_state = FindFinalState(current_settings,found);
                if(found){
                    ///if applicable, the found state will be reutilized in the next iteration
                    state = found_state;
                    Log() << "Previous state recovered, skipping\n";
                    ///if the result was found, we move on to the next iteration of the loop
                    continue;
                }
                else{
                    Log() << "No previous state recovered, going on\n";
                }
            }


            ///if the present iteration is the first one or we explicitly DON'T want to reutilize
            ///the final state of the previous iteration as the initial state for this one,
            ///we can do it here
            if(i==0 || settings.scanning.pass_on==false){
                ///this will be the present simulation, with the modified settings
                PRE79Simulation simulation(current_settings);
                ///this currently has no use. here we start a server which runs in a separate thread,
                ///which could possibly export information about the properties of the system at runtime
                if(settings.output.start_service){
                    RuntimePropertiesServer rt(current_settings,simulation);
                    ///an example for the usage of boost::thread
                    boost::thread thread1(rt);
                }
                
                ///"flush" the log from the initialization of the simulation to the present log
                Log() << simulation.GetLog();
                ///output will be written to the standard output
                simulation.SetStream(&Log());
                ///we have two options for the actual simulation, either entirely single-threaded
                ///or single-threaded thermalization followed by multi-threaded production
                ///both of these functions already contain thermalization, so we don't have to do it here
                if(settings.scanning.threaded_production)
                    ///run with parallel production (thermalization is still single-threaded)
                    state = simulation.RunParallel();
                else
                    ///run ordinary single-threaded simulation (both thermalization and production)
                    state = simulation.Run();
                    
            }
            ///here we will pass the final state of the previous iteration as an initial state for the present one
            else{
                ///this will be the present simulation, with the modified settings AND the "state" Lattice object as the initial state
                PRE79Simulation simulation(current_settings,state);
                ///this currently has no use. here we start a server which runs in a separate thread,
                ///which could possibly export information about the properties of the system at runtime
                if(settings.output.start_service){
                    RuntimePropertiesServer rt(current_settings,simulation);
                    ///an example for the usage of boost::thread
                    boost::thread thread1(rt);
                }
                
                ///"flush" the log from the initialization of the simulation to the present log
                Log() << simulation.GetLog();
                ///output will be written to the standard output
                simulation.SetStream(&Log());
                ///we have two options for the actual simulation, either entirely single-threaded
                ///or single-threaded thermalization followed by multi-threaded production
                ///both of these functions already contain thermalization, so we don't have to do it here
                if(settings.scanning.threaded_production)
                    ///run with parallel production (thermalization is still single-threaded)
                    state = simulation.RunParallel();
                else
                    ///run ordinary single-threaded simulation (both thermalization and production)
                    state = simulation.Run();
                    
            }
            ///between simulations, the remaining time is calculated
            pt::time_duration t1scan = pt::second_clock::local_time()-start_t;
            std::cout << "Total " << pt::to_simple_string(t1scan*(nscans-i-1)) << " remaining\n";
            start_t = pt::second_clock::local_time();
        }
        Log() << "Scanning finished\n";
    }
    /**
     * Naive parallelization, equal to executing many single-threaded simulations.
     * There is no obvious advantage to this method.
     */
    void RunParallel(){
        Log() << "Parallel OpenMP version\n";
        if(settings.scanning.separate_values)
            Log() << "Scanning values: " << settings.scanning.values << std::endl;
        else
            Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl;
        Log() << "Trying to set the number of concurrent simulations to " << settings.openmp.number_of_threads << std::endl;

        ///dynamic scheduling means that we don't know the actual number of threads prior to execution
        ///this isn't recommended for most purposes
        ///there seems to be some mechanism in OpenMP which chooses the appropriate number of threads
        ///but it doesn't seem to be very efficient
        if(settings.openmp.dynamic)
            omp_set_dynamic(1);
        else
            omp_set_dynamic(0);
        
        ///this function allows us to explicitly set number of threads
        omp_set_num_threads(settings.openmp.number_of_threads);
        ///automatic parallelization of the for loop, it runs nscans separate threads (not necessarily all threads at one time)
        ///schedule(runtime) means that we allow a runtime variable to change the scheduling, but
        ///actually we do it in main() through a small hack
        ///shared(rng2) is a common random number generator, which is used to initialize the thread-local random01
        ///this applies to the Mersenne-Twister case, but since 11.2011 it has been replaced by MWC
        ///and the directive here was left for compatibility
        ///private(random01) means that each thread has its own pseudorandom generator, initialized appropriately
        ///within the class (such as #rangen_mwc)
        #pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
        for(int i=0;i<nscans;i++){              ///run over all values of the parameters

            ///this is the actual present value of the parameter
            double value = 0.0;
            ///we set it to the appropriate value depending on the type of the set we are scanning over
            if(settings.scanning.separate_values)
                value = settings.scanning.values[i];
            else
                value = start + double(i)*delta;
            
            ///since we are multithreading, we need to make sure that non-thread-safe operations, such as working with streams, are executed one at a time
            ///pragma omp critical ensures that the following line will be executed by one thread at a time only
            #pragma omp critical
            Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ", Value: " << value << std::endl;

            ///we keep a copy of "localized" settings, where the actual present values are
            Settings current_settings(settings);

            ///several parameter names are supported
            ///we alter the "current settings" with accordance to the specified scanning variable name
            if(variable=="hamiltonian.tau")
                current_settings.hamiltonian.tau=value;
            if(variable=="hamiltonian.lambda")
                current_settings.hamiltonian.lambda=value;
            if(variable=="hamiltonian.temperature")
                current_settings.hamiltonian.temperature=value;
	    if(variable=="hamiltonian.h")
		current_settings.hamiltonian.h=value;

            ///here we search the database for an existing result for the present value
            if(current_settings.scanning.continue_if_results_exist){
                #pragma omp critical
                Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": Searching database for results of previous simulation\n";
                bool found=false;
                ///this function searches for the saved results and alters the "found" variable to true if it succeeds
                FindFinalProperties(current_settings,found);
                if(found){
                    #pragma omp critical
                    Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": Previous results recovered, skipping\n";
                    ///if the result was found, we move on to the next iteration of the loop
                    continue;
                }
                else{
                    #pragma omp critical
                    Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": No previous results recovered, going on\n";
                }
            }
            
            ///the actual simulation for the present thread
            PRE79Simulation thread_sim(current_settings);
            
            ///flush the log of the initialization procedure to the main log
            #pragma omp critical
            Log() << thread_sim.GetLog();
            ///set the output to standard output
            thread_sim.SetStream(&Log());
            ///run the simulation
            thread_sim.Run();
            
        }
        Log() << "Scanning finished\n";
    }
    /**
     * Two-fold parallelization. The thermalizations are run first in parallel 
     * until all the states are thermalized. Then we run parallel productions
     * on the thermalized states one by one, one at a time. This technique ensures
     * that the available processors are used to the maximum at all times.
     * 
     * If the number of parameter values is not divisible by the number of threads,
     * there may be some time when only the number of processors equal to the 
     * remainder of this division is used.
     */
    void RunParallelParallel() {
        Log() << "Parallel thermalization with parallel productions\n";
        if(settings.scanning.separate_values)
            Log() << "Scanning values: " << settings.scanning.values << std::endl;
        else
            Log() << "Scanning " << variable << " from " << start << " to " << end << " with interval " << delta << std::endl;
        Log() << "Trying to set the number of concurrent simulations to " << settings.openmp.number_of_threads << std::endl;

        ///dynamic scheduling means that we don't know the actual number of threads prior to execution
        ///this isn't recommended for most purposes
        ///there seems to be some mechanism in OpenMP which chooses the appropriate number of threads
        ///but it doesn't seem to be very efficient
        if(settings.openmp.dynamic)
            omp_set_dynamic(1);
        else
            omp_set_dynamic(0);
        ///explicitly set the number of threads
        omp_set_num_threads(settings.openmp.number_of_threads);

        ///the number of scans can be localized
        int nscans = PRE79Scanning::nscans;
        if(settings.scanning.separate_values)
            nscans = settings.scanning.values.size();

        ///we proceed in chunks of thermalizations. each chunk consists of a number
        ///of thremalizations equal at most to the number of threads
        ///inside this chunk, we first thermalize the states and then proceed to productions
        ///we then move on to the next chunk, thermalize the states, produce and so on
        ///incr is the size of the chunk
        int incr=settings.openmp.number_of_threads;
        ///note that this loop is not parallelized
        for(int i=0; i<nscans; i+=incr){

            int chunk = i/incr+1;
            int chunks = nscans/incr+1;

            int current_nscans = incr;
            if(chunk*incr > nscans)
                current_nscans = nscans-(chunk-1)*incr;
            
            ///a vector which will hold the thermalized states
            std::vector<Lattice>    states(current_nscans);


            Log() << "Chunk " << chunk << "/" << chunks << ": starting\n";
            Log() << "Chunk " << chunk << "/" << chunks << ": that is " << current_nscans << " thermalizations/productions of " << nscans << " total\n";
            Log() << "Chunk " << chunk << "/" << chunks << ": thermalization\n";

            ///parallel thermalization
            #pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
            for(int t=0;t<current_nscans;t++){
                
                ///this will be the actual value of the scanning parameter
                double value = 0.0;
                ///we will set it in accordance to the specified method of scannning, preset list of values or iteration over an interval
                if(settings.scanning.separate_values)
                    value = settings.scanning.values[t];
                else
                    value = start + double(t)*delta;
                
                ///make sure one thread at a time executes operations that are non-thread-safe
                #pragma omp critical
                Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ", Value: " << value << std::endl;

                ///a "local" copy of the settings we will pass on to the actual simulation
                Settings current_settings(settings);

                
                ///based on the specified name of the parameter, we alter this parameter in the "local" settings copy
                if(variable=="hamiltonian.tau")
                    current_settings.hamiltonian.tau=value;
                if(variable=="hamiltonian.lambda")
                    current_settings.hamiltonian.lambda=value;
                if(variable=="hamiltonian.temperature")
                    current_settings.hamiltonian.temperature=value;
                
                ///we look for a thermalized state that can be utilized in the next step,
                ///so we can omit this thermalization and just use the found state
                bool found = false;
                int cycle = 0;
                
                ///find the state, if not, who cares, we will set states[t] later anyway
                states[t] = FindLastStateTemperatureTolerant(current_settings,found,cycle);

                if(found){
                    #pragma omp critical
                    Log() << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() << ", Thermalized state found for value " << value << std::endl;
                    ///abort this thermalization, a thermalized state has already been found
                    continue;
                }
                
                ///no thermalized state was found, proceed
                ///create a thermalization simulation with localized settings
                PRE79Simulation termo(current_settings);
                ///flush the initialization log to the main log
                #pragma omp critical
                Log() << termo.GetLog();
                ///set the output to the main log
                termo.SetStream(&Log());
                ///thermalize the state and put it in our vector
                states[t]=termo.Thermalize();
            }

            Log() << "Chunk " << chunk << "/" << chunks << ": finished thermalization\n";
            Log() << "Chunk " << chunk << "/" << chunks << ": starting production\n";

            ///iterate AGAIN over the same chunk, this time for production
            for(int t=0;t<current_nscans;t++){
                ///the actual value of the scanning parameter
                double value = 0.0;
                ///we set it according to the chosen method of scannning
                if(settings.scanning.separate_values)
                    value = settings.scanning.values[t];
                else
                    value = start + double(t)*delta;
                
                Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ", Value: " << value << std::endl;

                ///localized settings copy
                Settings current_settings(settings);

                
                ///set the appropriate value in the localized copy of the settings according to the specified name of the scanning parameter
                if(variable=="hamiltonian.tau")
                    current_settings.hamiltonian.tau=value;
                if(variable=="hamiltonian.lambda")
                    current_settings.hamiltonian.lambda=value;
                if(variable=="hamiltonian.temperature")
                    current_settings.hamiltonian.temperature=value;

                ///entirely skip thermalization, since we have thermalized states already
                current_settings.simulation.thermalization_cycles=0;
                current_settings.simulation.supplementary_thermalization_cycles=0;

                ///here we look for existing results for the present value and if they exist, we skip the present production
                if(current_settings.scanning.continue_if_results_exist){
                    Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": Searching database for results of previous simulation\n";
                    bool found=false;
                    ///this function searches the database for existing results and alters the "found" variable to true if it succeeds
                    FindFinalProperties(current_settings,found);
                    if(found){
                        Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": Previous results recovered, skipping\n";
                        ///if we find results in the database, we omit the present calculation and skip to the next parameter value
                        continue;
                    }
                    else{
                        Log() << "Thread: "<< omp_get_thread_num() << "/" << omp_get_num_threads() << ": No previous results recovered, going on\n";
                    }
                }
                
                ///production for the local set of parameters, initialized with the state states[t]
                PRE79Simulation production(current_settings,states[t]);
                ///flush the initialization log to the main log
                Log() << production.GetLog();
                ///redirect the log to the main log
                production.SetStream(&Log());
                ///execute a parallel production for the present parameter values
                production.RunParallel();
            }
            Log() << "Chunk " << chunk << "/" << chunks << ": finished\n";
        }
        Log() << "Scanning finished\n";
    }

};


#endif	/* _PRE79SCANNING_H */


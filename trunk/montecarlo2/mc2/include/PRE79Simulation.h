/**
 * @file PRE79Simulation.h
 * 
 * Simulation classes
 * 
 */

#ifndef _PRE79SIMULATION_H
#define	_PRE79SIMULATION_H

#include "Simulation.h"
#include "serializer.h"
#include "ILoggable.h"
#include "SimulationDBFind.h"
#include "AutoCorrelationTimeCalculator.h"
#include <omp.h>

/**
 * A redundant class for a single production when a parallel production is performed.
 */
class PRE79Production:public ILoggable {
    shared_ptr<Lattice>                     lattice;            ///<pointer to the lattice object
    shared_ptr<PRE79StandardHamiltonian>    H;                  ///<pointer to the Hamiltonian object
    shared_ptr<LatticeSimulation>           simulation;         ///<pointer to the LatticeSimulation object
    shared_ptr<PRE79StandardProperties>     prop;               ///<pointer to the properties calculator and accumulator
    shared_ptr<Metropolis>                  metro;              ///<pointer to the Metropolis method
    Settings  &   settings;     ///<local copy of settings
    long nprod;                 ///<total number of production cycles
    long ncycles;               ///<number of cycles considered for measurement
public:
    /**
     * Constructor
     * 
     * @param set       reference to global Settings object
     * @param _nprod    number of production cycles for the entire production
     * @param start     initial state for production
     */
    PRE79Production(Settings & set,long _nprod, const Lattice & start):
    nprod(_nprod),
    settings(set)
    {
        ncycles = nprod/settings.simulation.measure_frequency;

        ///create the needed objects
        lattice = shared_ptr<Lattice>(new Lattice(start));
        H = shared_ptr<PRE79StandardHamiltonian>(new PRE79StandardHamiltonian(settings.hamiltonian.temperature, settings.hamiltonian.lambda, settings.hamiltonian.tau,settings.hamiltonian.h));
        metro = shared_ptr<Metropolis>(new Metropolis(settings,H,0.065));
        prop = shared_ptr<PRE79StandardProperties>(new PRE79StandardProperties(lattice,ncycles));
        simulation = shared_ptr<LatticeSimulation>(new LatticeSimulation(H,lattice,metro,nprod,0));
        //SetStream(&std::cout);
    }

    /**
     * Run a single production in the set.
     */
    void Run() {
        
#pragma omp critical
        Log() << "Production with freq " << settings.simulation.measure_frequency << std::endl ;

        ///save the time for calculation of remaining time
        pt::ptime start_t = pt::second_clock::local_time();
        ///the remaining time will be reported 5 times through the production
        int remaining_interval = simulation->GetNCycles()/5;
	if(settings.simulation.calculate_time)
        	Log() << "Remaining time will be reported every " << remaining_interval << " cycles\n";
        
        ///an instance of the autocorrelation time calculator
        shared_ptr<AutoCorrelationTimeCalculator> ac(new AutoCorrelationTimeCalculator(lattice,settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length));
        
        ///k will count the cycles
        long k=0;
        ///this loop will terminate if the underlying #LatticeSimulation Iterate() 
        ///funcion returns false, which will happen when all the cycles are finished
        ///the Iterate() function does the actual lattice sweeps
        while(simulation->Iterate()){
            ///update the autocorrelation time calculator
            ac->Update();
            
            ///do measurements here
            if(k%settings.simulation.measure_frequency==0){
                ///measurements are done via the Update() function of the properties object
                prop->Update(k,H,ac);
            }
            
            ///here the progress as percentage is calculated and reported
            if(k%1000==0 && settings.output.report_progress){
                Log() << "E = " << prop->EnergyEvolution()[k/1001] << std::endl;
                Log() << "Progress: " << (double(k)/double(simulation->GetNCycles()))*100.0 << "%\n";
            }
            
            
            ///runtime adjustment of the random-walk radius
            // Swendson (2011) the interval shouldn't be less than (total number of MC cycles)^(1/2)
            if(k%settings.simulation.radius_adjustment_frequency==0)
                metro->AdjustRadius(lattice);
            
            ///measure the acceptance rate if needed
            if(settings.simulation.measure_acceptance && k%settings.simulation.measure_acceptance_frequency==0){
                //Log() << "Acceptance rate: " << metro->MeasureAccepted(lattice)*100.0 << "%\n";
                Log() << "Acceptance rate: " << simulation->GetAcceptance()*100.0 << "%\n";
                Log() << "Mean acceptance rate: " << simulation->GetMeanAcceptance()*100.0 << "%\n";
            }

            ///calculate the remaining time of the present production
            if((k+1)%remaining_interval==0 && settings.simulation.calculate_time){
                pt::time_duration run1k = pt::second_clock::local_time() - start_t;
                int total1kruns = (simulation->GetNCycles()-k-1)/remaining_interval;
                std::cout << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() <<  ": "<< pt::to_simple_string(run1k*total1kruns) << " remaining\n";
                start_t = pt::second_clock::local_time();
            }
            k++;
        }
    }

    ///Access the pointer to the lattice object
    shared_ptr<Lattice> GetLattice(){
        return lattice;
    }
    ///Access the pointer to the properties object
    shared_ptr<PRE79StandardProperties> GetProperties(){
        return prop;
    }
    ///Access the pointer to the LatticeSimulation object (read only)
    const shared_ptr<LatticeSimulation> GetSimulation() const {
        return simulation;
    }
    /**
     * Redefinition of the Log stream, which adds thread information to the output.
     * The implementation at #ILoggable doesn't have this.
     */
    virtual std::ostream & Log(){
        ILoggable::Log() << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() << ": ";
	return ILoggable::Log();
    }
    ///This is actually unnecessary?
    void SetStream(std::ostream * os) {
	ILoggable::SetStream(os);
    }

};

/**
 * This class takes care of the entire process of preparing the simulation, thermalizing the initial states
 * and the actual production phase. Several mechanisms are implemented, such as
 * searching for thermalized states.
 */
class PRE79Simulation:public ILoggable {
    shared_ptr<Lattice>                     lattice;            ///<pointer to the lattice instance
    shared_ptr<PRE79StandardHamiltonian>    H;                  ///<pointer to the hamiltonian object
    shared_ptr<LatticeSimulation>           thermalization;     ///<underlying thermalization phase
    shared_ptr<LatticeSimulation>           simulation;         ///<underlying production phase
    shared_ptr<PRE79StandardProperties>     prop;               ///<pointer to production properties calculator
    shared_ptr<Metropolis>                  metro;              ///<pointer to the Metropolis method
    std::vector<shared_ptr<PRE79Production> >    productions;   ///<a vector of productions if parallel production is enabled
    shared_ptr<PRE79StandardProperties>     thermalprops;       ///<pointer to thermalization properties calculator
    Settings  &   settings;     ///<local settings reference
    SimulationDB        database;       ///<simulation database interface
    bool    restored;      ///<should we search for a thermalized state to reuse as an initial state?

    /**
     * Initialization
     */
    void Init(){
        Log() << "Creating Hamiltonian\n";
        H = shared_ptr<PRE79StandardHamiltonian>(new PRE79StandardHamiltonian(settings.hamiltonian.temperature, settings.hamiltonian.lambda, settings.hamiltonian.tau,settings.hamiltonian.h));
        Log() << "Creating Metropolis\n";
        metro = shared_ptr<Metropolis>(new Metropolis(settings,H,0.065));
        
        if(!restored){

            //--- szukamy ewentualnego zapisanego stermalizowanego stanu
            if(settings.simulation.find_thermalized) {
                Log() << "Searching database for thermalized state\n";
                bool found=false;

                //ostatni zapisany stan lepiej nadaje się jako stermalizowany
                int cycle=0; // nieużywane

                //--- dopuszczamy załadowanie stanu z zadaną tolerancją temperatury, ale wtedy musimy dotermalizować
                Lattice found_state;
                if(settings.simulation.find_thermalized_temperature_tolerance!=0.0) {
                    Log() << "Temperature tolerance enabled\n";
                    found_state = FindLastStateTemperatureTolerant(settings,found,cycle,settings.simulation.find_thermalized_temperature_tolerance);
                }
                else {
                    found_state = FindLastState(settings,found,cycle);
                }
                //---
                //--- dopuszczamy ewentualny stan z polem o bliskiej wartości, ale musimy dotermalizować
                if(settings.simulation.find_thermalized_h_tolerance!=0.0) {
                    Log() << "Field tolerance enabled\n";
                    found_state = FindLastStateFieldTolerant(settings,found,cycle,settings.simulation.find_thermalized_h_tolerance);
                }
                else {
                    found_state = FindLastState(settings,found,cycle);
                }
                //---

                if(found){
                    lattice = shared_ptr<Lattice>(new Lattice(found_state));
                    
                    Log() << "Creating Supplementary Thermalization for recovered state\n";
                    thermalization = shared_ptr<LatticeSimulation>(new LatticeSimulation(H,lattice,metro,settings.simulation.supplementary_thermalization_cycles));
                }
                else {
                    Log() << "No thermalized state found, creating Thermalization\n";
                    thermalization = shared_ptr<LatticeSimulation>(new LatticeSimulation(H,lattice,metro,settings.simulation.thermalization_cycles));
                }
            }
            //---
            else {
                Log() << "No thermalized state found, creating Thermalization\n";
                thermalization = shared_ptr<LatticeSimulation>(new LatticeSimulation(H,lattice,metro,settings.simulation.thermalization_cycles));
            }
        }
        else{
            Log() << "Previous final state passed on as initial, creating Supplementary Thermalization\n";
            thermalization = shared_ptr<LatticeSimulation>(new LatticeSimulation(H,lattice,metro,settings.simulation.supplementary_thermalization_cycles));
        }
        //--- szukamy stanu, od którego możemy kontynuować symulację
        //TODO: stan się ładuje, ale jeszcze trzeba policzyć albo wczytać brakujące Properties...
        int found_cycle = 0;
        if(settings.simulation.pick_up_aborted){
            Log() << "Searching database for aborted simulation\n";
            
            bool found = false;
            Lattice found_state = FindLastState(settings,found,found_cycle);
            if(found){
                
                lattice = shared_ptr<Lattice>(new Lattice(found_state));
                Log() << "Recovering from aborted simulation at " << found_cycle << " cycles\n";
            } else {
                Log() << "No aborted simulation found\n";
            }
        }
        //---
        int cycle_advantage = settings.simulation.production_cycles - found_cycle;

        Log() << "Creating Production\n";
        simulation = shared_ptr<LatticeSimulation>(new LatticeSimulation(H,lattice,metro,cycle_advantage,found_cycle));
        Log() << "Creating Properties\n";
        prop = shared_ptr<PRE79StandardProperties>(new PRE79StandardProperties(lattice,settings.simulation.production_cycles/settings.simulation.measure_frequency));
        thermalprops = shared_ptr<PRE79StandardProperties>(new PRE79StandardProperties(lattice,thermalization->GetNCycles()/10));

    }

public:
    PRE79Simulation(Settings & set):
    settings(set),
    database(set),
            thermalprops(shared_ptr<PRE79StandardProperties>()),thermalization(shared_ptr<LatticeSimulation>()),lattice(shared_ptr<Lattice>())
    {
        if(set.initial.isotropic)
            lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::Isotropic));
        else
        if(set.initial.biaxial){
            if(set.initial.righthanded)
                lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::BiaxialRighthanded));
            else
                lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::Biaxial));
        }
        else
        if(set.initial.biaxial_alt){
            if(set.initial.righthanded)
                lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::BiaxialRighthandedAlt));
            else
                lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::BiaxialAlt));
        }
        else
        if(set.initial.righthanded)
            lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::IsotropicRighthanded));
        else
            lattice = shared_ptr<Lattice>(new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::Isotropic));
        restored=false;
        Init();
    }
    PRE79Simulation(Settings & set,const Lattice & saved):
    settings(set),database(set)
    {
        //nie kopiujemy wskaźnika, kopiujemy obiekt
        //żeby nie zmieniać źródłowego obiektu
        restored=true;
        lattice=shared_ptr<Lattice>(new Lattice(saved));
        Init();
    }
    const pt::time_duration DurationOf1000Cycles() const {
        LatticeSimulation test(H,lattice,metro,1000);
        PRE79StandardProperties testprop(lattice,1000/settings.simulation.measure_frequency);
        pt::ptime start_t = pt::second_clock::local_time();

        int k=0;
        while(test.Iterate()){
            if(k%settings.simulation.measure_frequency==0)
                testprop.Update(k,H);
            k++;
        }
        pt::time_duration run1k = pt::second_clock::local_time() - start_t;
        return run1k;
    }
    const pt::time_duration ExpectedSimulationTime() const {
        
        int total1kruns = (thermalization->GetNCycles()+simulation->GetNCycles())/1000 ;
        return DurationOf1000Cycles()*total1kruns;
    }

    const Lattice & Thermalize() {
        shared_ptr<AutoCorrelationTimeCalculator> ac(new AutoCorrelationTimeCalculator(lattice,settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length));
        Log() << "Thermalization cycles: " << thermalization->GetNCycles() << std::endl;
        int tcycle=0;
        while(thermalization->Iterate()){
            ac->Update();
            if(tcycle%10==0)
                thermalprops->Update(tcycle,H,ac);
            if(tcycle%1000==0 && settings.output.report_progress){
                Log() << "E = " << thermalprops->EnergyEvolution()[tcycle/1001] << std::endl;
                Log() << "Progress: " << (double(tcycle)/double(thermalization->GetNCycles()))*100.0 << "%\n";
            }
            //--- poprawa promienia błądzenia przypadkowego <-- czyżby źródło błędów???
            if(tcycle%settings.simulation.radius_adjustment_frequency==0)
                metro->AdjustRadius(lattice);
            //---
            if(settings.simulation.measure_acceptance && tcycle%settings.simulation.measure_acceptance_frequency==0){
                //Log() << "Acceptance rate: " << metro->MeasureAccepted(lattice)*100.0 << "%\n";
                Log() << "Acceptance rate: " << thermalization->GetAcceptance()*100.0 << "%\n";
                Log() << "Mean acceptance rate: " << thermalization->GetMeanAcceptance()*100.0 << "%\n";
            }
            tcycle++;
        }
        //-- zapisywanie historii termalizacji
        
        if(settings.output.save_thermalization_properties){
            Log() << "Saving thermalization history\n";
            database.StoreThermalizationHistory(settings,*thermalprops);
        }

        return *lattice;
    }

    ///jednowątkowa termalizacja i wielowątkowa produkcja
    Lattice RunParallel(){
        //--- termalizacja
        Log() << "Adjusting radius\n";
        metro->AdjustRadius(lattice);
        if(settings.simulation.calculate_time){
            Log() << "Calculating expected duration of simulation\n";
            pt::time_duration   expected = DurationOf1000Cycles()*(thermalization->GetNCycles()+settings.simulation.production_cycles/settings.openmp.number_of_threads)/1000;
            Log() << "Expected time of simulation: " << pt::to_simple_string(expected) << std::endl;
        }
        if(settings.simulation.thermalization_cycles>0 && settings.simulation.supplementary_thermalization_cycles>0){
            Log() << "Thermalization\n";
            Thermalize();
        }
        //---
        
        //--- tworzymy fragmentaryczne symulacje
        
        //tutaj nie może być wątpliwości ile jest wątków
        omp_set_dynamic(0);
        omp_set_num_threads(settings.openmp.number_of_threads);
        
        //PRE79Production prototype(settings,settings.simulation.production_cycles/settings.openmp.number_of_threads,*lattice);
        //(settings.openmp.number_of_threads,prototype);
        for(int i=0;i<settings.openmp.number_of_threads;i++)
            productions.push_back(shared_ptr<PRE79Production>(new PRE79Production(settings,settings.simulation.production_cycles/settings.openmp.number_of_threads,*lattice)));

        Log() << "Starting " << settings.openmp.number_of_threads << " productions\n";
        #pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
        for(int i=0;i<productions.size();i++){
	    productions[i]->SetStream(&Log());
            productions[i]->Run();
        }
        shared_ptr<PRE79StandardProperties> generalprop = productions[0]->GetProperties();
        for(int i=1;i<productions.size();i++){
            generalprop->Append(productions[i]->GetProperties());
        }
        generalprop->CalculateSpecificHeat();

        if(settings.output.save_properties_evolution) {
            Log() << "Saving Properties Evolution\n";
            database.StorePropertiesEvolution(settings,*generalprop);
        }
        if(settings.output.save_final_configuration){
            Log() << "Saving Final Lattice\n";
            database.StoreFinalLattice(settings,*productions[0]->GetLattice());
        }
        if(settings.output.save_final_properties){
            Log() << "Saving Final Properties\n";
            PRE79MeanProperties pp(generalprop,H);
            database.StoreFinalProperties(settings,pp);
        }

        Log() << "Done\n";
        Log() << "Mean EPM: " << generalprop->TemporalMeanEnergyPerMolecule().Print() << std::endl;
        Log() << "Specific Heat: " << generalprop->SpecificHeat().Print() << std::endl;

        //return the state with the lowest energy
        vect weights(0.0,productions.size());
        for(int i=0;i<weights.size();i++){
            weights[i]=productions[i]->GetLattice()->GetMeanEPM();
        }
        int best = MinimumIndex(weights);

        Lattice ret = *productions[best]->GetLattice();
        //--

        productions.clear();

        return ret;
    }


    ///Symulacja. Zwraca końcowy stan sieci.
    const Lattice & Run(){
        Log() << "Production cycles: " << simulation->GetNCycles() << std::endl;
        Log() << "Adjusting radius\n";
        metro->AdjustRadius(lattice);
        if(settings.simulation.calculate_time){
            Log() << "Calculating expected duration of simulation\n";
            pt::time_duration   expected=ExpectedSimulationTime();
            Log() << "Expected time of simulation: " << pt::to_simple_string(expected) << std::endl;
        }
        if(settings.simulation.thermalization_cycles>0 && settings.simulation.supplementary_thermalization_cycles>0){
            Log() << "Thermalization\n";
            Thermalize();
        }

        //--- zapisywanie stanów pośrednich
        int intermediate_frequency = simulation->GetNCycles()/settings.output.intermediate_states;
        //---

        #pragma omp critical
        Log() << "Production with freq " << settings.simulation.measure_frequency << std::endl ;
        long k=0;

        pt::ptime start_t = pt::second_clock::local_time();
        int remaining_interval = simulation->GetNCycles()/20+1;
	if(settings.simulation.calculate_time)
	    Log() << "Remaining time will be reported every " << remaining_interval << " cycles\n";
        
        shared_ptr<AutoCorrelationTimeCalculator> ac(new AutoCorrelationTimeCalculator(lattice,settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length));
        
        while(simulation->Iterate()){
            ac->Update();
            //--- pomiary
            if(k%settings.simulation.measure_frequency==0){
                prop->Update(k,H,ac);
                if(settings.output.save_configuration_evolution){
                    database.StoreLattice(settings,*lattice,k);
                }
            }
            //---
            
            //--- sprawozdanie z postępu symulacji
            if(k%1000==0 && settings.output.report_progress){
                Log() << "E = " << prop->EnergyEvolution()[k/1001] << std::endl;
                Log() << "Progress: " << (double(k)/double(simulation->GetNCycles()))*100.0 << "%\n";
            }
            //---

            //--- zapis stanów pośrednich
            if(k%intermediate_frequency==0){
                if(settings.output.save_intermediate_states){
                    database.StoreLattice(settings,*lattice,k);
                    PRE79MeanProperties pp(prop,H);
                    database.StoreProperties(settings,pp,k);
                }

            }

            //--- poprawa promienia błądzenia przypadkowego <-- czyżby źródło błędów???
            // Swendson (2011) częstość zmiany kroku Monte Carlo nie powinna być niższa, niż (liczba kroków MC)^(1/2)
            if(k%settings.simulation.radius_adjustment_frequency==0)
                metro->AdjustRadius(lattice);
            //---

            //--- pomiar poziomu akceptacji
            if(settings.simulation.measure_acceptance && k%settings.simulation.measure_acceptance_frequency==0){
                //Log() << "Acceptance rate: " << metro->MeasureAccepted(lattice)*100.0 << "%\n";
                Log() << "Acceptance rate: " << simulation->GetAcceptance()*100.0 << "%\n";
                Log() << "Mean acceptance rate: " << simulation->GetMeanAcceptance()*100.0 << "%\n";
            }
            //---
            
            if((k+1)%remaining_interval==0 && settings.simulation.calculate_time){
                pt::time_duration run1k = pt::second_clock::local_time() - start_t;
                int total1kruns = (simulation->GetNCycles()-k-1)/remaining_interval;
                std::cout << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() <<  ": "<< pt::to_simple_string(run1k*total1kruns) << " remaining\n";
                start_t = pt::second_clock::local_time();
            }
            k++;
        }

        
        if(settings.output.save_properties_evolution) {
            Log() << "Saving Properties Evolution\n";
            database.StorePropertiesEvolution(settings,*prop);
        }
        if(settings.output.save_final_configuration){
            Log() << "Saving Final Lattice\n";
            database.StoreFinalLattice(settings,*lattice);
        }
        if(settings.output.save_final_properties){
            Log() << "Saving Final Properties\n";
            PRE79MeanProperties pp(prop,H);
            database.StoreFinalProperties(settings,pp);
        }
        
        Log() << "Done\n";
        Log() << "Mean EPM: " << prop->TemporalMeanEnergyPerMolecule().Print() << std::endl;
        Log() << "Specific Heat: " << prop->SpecificHeat().Print() << std::endl;
        return *lattice;

    }
    virtual std::ostream & Log(){
        ILoggable::Log() << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() << ": ";
	return ILoggable::Log();
    }
    /*
    
    ///for parallel tempering
    double SwapTemperature(const double & t) {
        double oldT = H->GetTemperature();
        double oldL = H->GetLambda();
        double oldTau = H->GetTau();
        double oldH = H->GetH();
        
        delete H;
        H = new PRE79StandardHamiltonian(t,oldL,oldTau,oldH);
        return oldT;
    }
    */

    shared_ptr<Lattice> GetLattice(){
        return lattice;
    }
    shared_ptr<PRE79StandardProperties> GetProperties(){
        return prop;
    }
    shared_ptr<PRE79StandardProperties> GetThermalizationProperties(){
        return thermalprops;
    }
    int GetNProductions(){
        return productions.size();
    }
    shared_ptr<PRE79Production> GetProduction(const int & i) {
        return productions[i];
    }
    shared_ptr<LatticeSimulation> GetSimulation() const {
        return simulation;
    }
    shared_ptr<LatticeSimulation> GetThermalization() const {
        return thermalization;
    }

};


#endif	/* _PRE79SIMULATION_H */


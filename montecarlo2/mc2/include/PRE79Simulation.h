/* 
 * File:   PRE79Simulation.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 19:47
 */

#ifndef _PRE79SIMULATION_H
#define	_PRE79SIMULATION_H

#include "Simulation.h"
#include "serializer.h"
#include "ILoggable.h"
#include "SimulationDBFind.h"
#include <omp.h>

///Podstawowa symulacja (produkcja) do produkcji równoległej
class PRE79Production:public ILoggable {
    Lattice                     *lattice;
    PRE79StandardHamiltonian    *H;
    LatticeSimulation           *simulation;
    PRE79StandardProperties     *prop;
    Metropolis                  *metro;
    const Settings  &   settings;
    long nprod;
    long ncycles;
public:
    PRE79Production(const Settings & set,long _nprod, const Lattice & start):
    nprod(_nprod),
    settings(set)
    {
        ncycles = nprod/settings.simulation.measure_frequency;

        lattice = new Lattice(start);
        H = new PRE79StandardHamiltonian(settings.hamiltonian.temperature, settings.hamiltonian.lambda, settings.hamiltonian.tau,settings.hamiltonian.h);
        metro = new Metropolis(settings,H,0.065);
        prop = new PRE79StandardProperties(lattice,ncycles);
        simulation = new LatticeSimulation(H,lattice,metro,nprod,0);
        //SetStream(&std::cout);
    }

    void Run() {
        //std::cout << "Thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
        #pragma omp critical
        Log() << "Production with freq " << settings.simulation.measure_frequency << std::endl ;



        pt::ptime start_t = pt::second_clock::local_time();
        int remaining_interval = simulation->GetNCycles()/5;
	if(settings.simulation.calculate_time)
        	Log() << "Remaining time will be reported every " << remaining_interval << " cycles\n";

        long k=0;
        while(simulation->Iterate()){
            //--- pomiary
            if(k%settings.simulation.measure_frequency==0){
                prop->Update(k,H);
            }
            //---
            //--- poprawa promienia błądzenia przypadkowego <-- czyżby źródło błędów???
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
    }

    Lattice & GetLattice(){
        return *lattice;
    }
    PRE79StandardProperties & GetProperties(){
        return *prop;
    }
    const LatticeSimulation * GetSimulation() const {
        return simulation;
    }
    virtual std::ostream & Log(){
        ILoggable::Log() << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() << ": ";
	return ILoggable::Log();
    }
    void SetStream(std::ostream * os) {
	ILoggable::SetStream(os);
    }
    ~PRE79Production(){
        delete lattice;
        delete H;
        delete metro;
        delete prop;
        delete simulation;
    }
};

///Pojedyncza symulacja (bez skanowania)
class PRE79Simulation:public ILoggable {
    Lattice                     *lattice;
    PRE79StandardHamiltonian    *H;
    LatticeSimulation           *thermalization;
    LatticeSimulation           *simulation;
    PRE79StandardProperties     *prop;
    Metropolis                  *metro;
    std::vector<PRE79Production*>    productions;
    PRE79StandardProperties     *thermalprops;
    const Settings  &   settings;
    SimulationDB        database;
    bool    restored;      ///<zaczynamy z wczytanego stanu sieci

    void Init(){
        Log() << "Creating Hamiltonian\n";
        H = new PRE79StandardHamiltonian(settings.hamiltonian.temperature, settings.hamiltonian.lambda, settings.hamiltonian.tau,settings.hamiltonian.h);
        Log() << "Creating Metropolis\n";
        metro = new Metropolis(settings,H,0.065);
        //metro->SetStream(&std::cout);
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
                    delete lattice;
                    lattice = new Lattice(found_state);
                    
                    Log() << "Creating Supplementary Thermalization for recovered state\n";
                    thermalization = new LatticeSimulation(H,lattice,metro,settings.simulation.supplementary_thermalization_cycles);
                }
                else {
                    Log() << "No thermalized state found, creating Thermalization\n";
                    thermalization = new LatticeSimulation(H,lattice,metro,settings.simulation.thermalization_cycles);
                }
            }
            //---
            else {
                Log() << "No thermalized state found, creating Thermalization\n";
                thermalization = new LatticeSimulation(H,lattice,metro,settings.simulation.thermalization_cycles);
            }
        }
        else{
            Log() << "Previous final state passed on as initial, creating Supplementary Thermalization\n";
            thermalization = new LatticeSimulation(H,lattice,metro,settings.simulation.supplementary_thermalization_cycles);
        }
        //--- szukamy stanu, od którego możemy kontynuować symulację
        //TODO: stan się ładuje, ale jeszcze trzeba policzyć albo wczytać brakujące Properties...
        int found_cycle = 0;
        if(settings.simulation.pick_up_aborted){
            Log() << "Searching database for aborted simulation\n";
            
            bool found = false;
            Lattice found_state = FindLastState(settings,found,found_cycle);
            if(found){
                delete lattice;
                lattice = new Lattice(found_state);
                Log() << "Recovering from aborted simulation at " << found_cycle << " cycles\n";
            } else {
                Log() << "No aborted simulation found\n";
            }
        }
        //---
        int cycle_advantage = settings.simulation.production_cycles - found_cycle;

        Log() << "Creating Production\n";
        simulation = new LatticeSimulation(H,lattice,metro,cycle_advantage,found_cycle);
        Log() << "Creating Properties\n";
        prop = new PRE79StandardProperties(lattice,settings.simulation.production_cycles/settings.simulation.measure_frequency);
        thermalprops = new PRE79StandardProperties(lattice,thermalization->GetNCycles()/100);

    }

public:
    PRE79Simulation(const Settings & set):
    settings(set),
    database(set),
            thermalprops(NULL),thermalization(NULL),lattice(NULL)
    {
        if(set.initial.isotropic)
            lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::Isotropic);
        else
        if(set.initial.biaxial){
            if(set.initial.righthanded)
                lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::BiaxialRighthanded);
            else
                lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::Biaxial);
        }
        else
        if(set.initial.biaxial_alt){
            if(set.initial.righthanded)
                lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::BiaxialRighthandedAlt);
            else
                lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::BiaxialAlt);
        }
        else
        if(set.initial.righthanded)
            lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::IsotropicRighthanded);
        else
            lattice = new Lattice(set.lattice.L,set.lattice.W,set.lattice.H,Lattice::Isotropic);
        restored=false;
        Init();
    }
    PRE79Simulation(const Settings & set,const Lattice & saved):
    settings(set),database(set)
    {
        //nie kopiujemy wskaźnika, kopiujemy obiekt
        //żeby nie zmieniać źródłowego obiektu
        restored=true;
        lattice=new Lattice(saved);
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
        Log() << "Thermalization cycles: " << thermalization->GetNCycles() << std::endl;
        int tcycle=0;
        while(thermalization->Iterate()){
            if(tcycle%100==0)
                thermalprops->Update(tcycle,H);
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
	Log() << "Saving thermalization history\n";
        database.StoreThermalizationHistory(settings,*thermalprops);

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
            productions.push_back(new PRE79Production(settings,settings.simulation.production_cycles/settings.openmp.number_of_threads,*lattice));

        Log() << "Starting " << settings.openmp.number_of_threads << " productions\n";
        #pragma omp parallel for schedule(runtime) shared(rng2) private(random01)
        for(int i=0;i<productions.size();i++){
	    productions[i]->SetStream(&Log());
            productions[i]->Run();
        }
        PRE79StandardProperties generalprop = productions[0]->GetProperties();
        for(int i=1;i<productions.size();i++){
            generalprop.Append(productions[i]->GetProperties());
        }
        generalprop.CalculateSpecificHeat();

        if(settings.output.save_properties_evolution) {
            Log() << "Saving Properties Evolution\n";
            database.StorePropertiesEvolution(settings,generalprop);
        }
        if(settings.output.save_final_configuration){
            Log() << "Saving Final Lattice\n";
            database.StoreFinalLattice(settings,productions[0]->GetLattice());
        }
        if(settings.output.save_final_properties){
            Log() << "Saving Final Properties\n";
            PRE79MeanProperties pp(generalprop,*H);
            database.StoreFinalProperties(settings,pp);
        }

        Log() << "Done\n";
        Log() << "Mean EPM: " << generalprop.TemporalMeanEnergyPerMolecule().Print() << std::endl;
        Log() << "Specific Heat: " << generalprop.SpecificHeat().Print() << std::endl;

        //return the state with the lowest energy
        vect weights(0.0,productions.size());
        for(int i=0;i<weights.size();i++){
            weights[i]=productions[i]->GetLattice().GetMeanEPM();
        }
        int best = MinimumIndex(weights);

        Lattice ret = productions[best]->GetLattice();
        //--

        //cleanup
        foreach(PRE79Production * prod,productions){
            delete prod;
        }

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
        while(simulation->Iterate()){
            //--- pomiary
            if(k%settings.simulation.measure_frequency==0){
                prop->Update(k,H);
                if(settings.output.save_configuration_evolution){
                    database.StoreLattice(settings,*lattice,k);
                }
            }
            //---

            //--- zapis stanów pośrednich
            if(k%intermediate_frequency==0){
                if(settings.output.save_intermediate_states){
                    database.StoreLattice(settings,*lattice,k);
                    PRE79MeanProperties pp(*prop,*H);
                    database.StoreProperties(settings,pp,k);
                }

            }

            //--- poprawa promienia błądzenia przypadkowego <-- czyżby źródło błędów???
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
            PRE79MeanProperties pp(*prop,*H);
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
    ~PRE79Simulation(){
        delete metro;
        delete lattice;
        delete H;
        delete prop;
        delete thermalization;
        delete simulation;
        delete thermalprops;
    }

    Lattice & GetLattice(){
        return *lattice;
    }
    PRE79StandardProperties & GetProperties(){
        return *prop;
    }
    PRE79StandardProperties & GetThermalizationProperties(){
        return *thermalprops;
    }
    int GetNProductions(){
        return productions.size();
    }
    PRE79Production & GetProduction(const int & i) {
        return *productions[i];
    }
    const LatticeSimulation * GetSimulation() const {
        return simulation;
    }
    const LatticeSimulation * GetThermalization() const {
        return thermalization;
    }

};


#endif	/* _PRE79SIMULATION_H */


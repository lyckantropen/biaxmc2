/* 
 * File:   PRE79Simulation.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 19:47
 */

#ifndef _PRE79SIMULATION_H
#define	_PRE79SIMULATION_H

#include "Simulation.h"
#include "Settings.h"
#include "serializer.h"
#include "SimulationDB.h"
#include "PRE79StandardProperties.h"
#include "ILoggable.h"
#include "SimulationDBFind.h"
#include <omp.h>

///Pojedyncza symulacja (bez skanowania)
class PRE79Simulation:public ILoggable {
    Lattice                     *lattice;
    PRE79StandardHamiltonian    *H;
    LatticeSimulation           *thermalization;
    LatticeSimulation           *simulation;
    PRE79StandardProperties     *prop;
    Metropolis                  *metro;
    const Settings  &   settings;
    SimulationDB        database;
    bool    restored;      ///<zaczynamy z wczytanego stanu sieci

    void Init(){
        Log() << "Creating Hamiltonian\n";
        H = new PRE79StandardHamiltonian(settings.hamiltonian.temperature, settings.hamiltonian.lambda, settings.hamiltonian.tau);
        Log() << "Creating Metropolis\n";
        metro = new Metropolis(H,0.065);
        if(!restored){
            //--- szukamy ewentualnego zapisanego stermalizowanego stanu
            if(settings.simulation.find_thermalized) {
                Log() << "Searching database for thermalized state\n";
                bool found=false;
                //Lattice found_state = FindState(settings,0,found);
                //ostatni zapisany stan lepiej nadaje się jako stermalizowany
                int cycle=0; // nieużywane
                Lattice found_state = FindLastState(settings,found,cycle);
                if(found){
                    delete lattice;
                    lattice = new Lattice(found_state);
                    //pusta termalizacja
                    thermalization = new LatticeSimulation(H,lattice,metro,0);
                    Log() << "Thermalized state found, picking up\n";
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
            Log() << "Creating Supplementary Thermalization\n";
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
    }

public:
    PRE79Simulation(const Settings & set):
    settings(set),
    database(set)
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

    ///Symulacja. Zwraca końcowy stan sieci.
    const Lattice * Run(){
        Log() << "Thermalization cycles: " << thermalization->GetNCycles() << std::endl;
        Log() << "Production cycles: " << simulation->GetNCycles() << std::endl;
        Log() << "Adjusting radius\n";
        metro->AdjustRadius(lattice);
        Log() << "Calculating expected duration of simulation\n";
        pt::time_duration   expected=ExpectedSimulationTime();
        Log() << "Expected time of simulation: " << pt::to_simple_string(expected) << std::endl;
        Log() << "Thermalization\n";
        while(thermalization->Iterate());

        //--- zapisywanie stanów pośrednich
        int intermediate_frequency = simulation->GetNCycles()/settings.output.intermediate_states;
        //---

        Log() << "Production with freq " << settings.simulation.measure_frequency << std::endl ;
        int k=0;

        pt::ptime start_t = pt::second_clock::local_time();
        int remaining_interval = expected.total_seconds()*100.0/std::sqrt(thermalization->GetNCycles()/1000+simulation->GetNCycles()/1000);
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
                if(settings.output.save_intermediate_states)
                    database.StoreLattice(settings,*lattice,k);
            }

            //--- poprawa promienia błądzenia przypadkowego
            if(k%settings.simulation.radius_adjustment_frequency==0)
                metro->AdjustRadius(lattice);
            //---
            
            if((k+1)%remaining_interval==0){
                pt::time_duration run1k = pt::second_clock::local_time() - start_t;
                int total1kruns = (simulation->GetNCycles()-k-1)/remaining_interval;
                std::cout << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() <<  ": "<< pt::to_simple_string(run1k*total1kruns) << " remaining\n";
                start_t = pt::second_clock::local_time();
            }
            k++;
        }

        
        if(settings.output.save_properties_evolution)
            database.StorePropertiesEvolution(settings,*prop);
        if(settings.output.save_final_configuration)
            database.StoreFinalLattice(settings,*lattice);
        if(settings.output.save_final_properties){
            PRE79MeanProperties pp(*prop,*H);
            database.StoreFinalProperties(settings,pp);
        }
        
        Log() << "Done\n";
        Log() << "Mean EPM: " << prop->TemporalMeanEnergyPerMolecule().Print() << std::endl;
        Log() << "Specific Heat: " << prop->SpecificHeat().Print() << std::endl;
        Log() << "Uniaxial From Correlation: " << prop->UniaxialOrderByCorrelation().Print() << std::endl;
        Log() << "Biaxial From Correlation: " << prop->BiaxialOrderByCorrelation().Print() << std::endl;
        Log() << "Tetra From Correlation: " << prop->TetrahedralOrderByCorrelation().Print() << std::endl;
        Log() << "Parity From Correlation: " << prop->ParityOrderByCorrelation().Print() << std::endl;
        return lattice;

    }
    virtual std::ostream & Log(){
        ILoggable::Log() << "Thread: " << omp_get_thread_num() << "/" << omp_get_num_threads() << ": ";
    }
    ~PRE79Simulation(){
        delete metro;
        delete lattice;
        delete H;
        delete prop;
        delete thermalization;
        delete simulation;
    }
};


#endif	/* _PRE79SIMULATION_H */


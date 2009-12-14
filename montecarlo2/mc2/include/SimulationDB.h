/* 
 * File:   SimulationDB.h
 * Author: karol
 *
 * Created on 4 grudzień 2009, 12:12
 */

#ifndef _SIMULATIONDB_H
#define	_SIMULATIONDB_H

#include "base.h"
#include "boost.h"
#include "Settings.h"
#include "PRE79StandardProperties.h"
#include "PRE79SpatialCorrelations.h"

class SimulationDB {
    boostbase::base db;
    std::string id;         ///<identyfikator projektu
    std::string user;       ///<użytkownik uzyskujący dostęp do bazy
    int ncycles;

public:
    const std::string id_kw;
    const std::string type_label;
    const std::string L_label;
    const std::string W_label;
    const std::string H_label;
    const std::string user_label;
    const std::string final_properties_kw;  ///<słowo kluczowe końcowych właściwości systemu
    const std::string current_lattice_kw;           ///<słowo kluczowe znaczące stan siatki 
    const std::string final_lattice_kw;     ///<słowo kluczowe znaczące ostateczny stan siatki
    const std::string properties_evolution_kw;        ///<słowo kluczowe znaczące całość ewolucji właściwości systemu
    const std::string current_properties_kw;
    const std::string cycle_label;             ///<słowo kluczowe oznaczające numer cyklu
    const std::string ncycles_label;           ///<słowo kluczowe oznaczające ilość wszystkich cykli w symulacji
    const std::string temperature_label;
    const std::string lambda_label;
    const std::string tau_label;
public:
    SimulationDB(const Settings & set):
    db(fs::path(set.sqlite.file),fs::path(set.sqlite.dir)),
            L_label("L"),
            W_label("W"),
            H_label("H"),
            user_label("user"),
            type_label("data_type"),
            id_kw("project_name"),
            final_properties_kw("final_properties"),
            current_lattice_kw("lattice"),
            final_lattice_kw("final_lattice"),
            properties_evolution_kw("properties_evolution"),
            current_properties_kw("properties"),
            cycle_label("production_cycle"),
            ncycles_label("total_production_cycles"),
            temperature_label("temperature"),
            lambda_label("lambda"),
            tau_label("tau")
    {
        id = set.project.name;
        user = std::getenv("USER");
        ncycles = set.simulation.production_cycles;
    }
    void StoreLattice(const Settings & settings,Lattice & lat,const int & cycle){
        db.store<Lattice>(lat,boostbase::where
          (id_kw,id)
          (user_label,user)
          (L_label,settings.lattice.L)
          (W_label,settings.lattice.W)
          (H_label,settings.lattice.H)
          (cycle_label,cycle)
          (ncycles_label,ncycles)
          (temperature_label,settings.hamiltonian.temperature)
          (lambda_label,settings.hamiltonian.lambda)
          (tau_label,settings.hamiltonian.tau)
          (type_label,current_lattice_kw)
        );
    }
    void StoreFinalLattice(const Settings & settings,Lattice & lat){
        db.store<Lattice>(lat,boostbase::where
          (id_kw,id)
          (user_label,user)
          (L_label,settings.lattice.L)
          (W_label,settings.lattice.W)
          (H_label,settings.lattice.H)
          (cycle_label,(ncycles-1))
          (ncycles_label,ncycles)
          (temperature_label,settings.hamiltonian.temperature)
          (lambda_label,settings.hamiltonian.lambda)
          (tau_label,settings.hamiltonian.tau)
          (type_label,final_lattice_kw)
        );
    }
    void StorePropertiesEvolution(const Settings & settings,PRE79StandardProperties & prop){
        db.store<PRE79StandardProperties>(prop,boostbase::where
          (id_kw,id)
          (user_label,user)
          (L_label,settings.lattice.L)
          (W_label,settings.lattice.W)
          (H_label,settings.lattice.H)
          (cycle_label,ncycles-1)
          (ncycles_label,ncycles)
          (temperature_label,settings.hamiltonian.temperature)
          (lambda_label,settings.hamiltonian.lambda)
          (tau_label,settings.hamiltonian.tau)
          (type_label,properties_evolution_kw)
        );
    }
    const std::stringstream & SqliteLog(){
        return db.log();
    }

    void StoreProperties(const Settings & settings,PRE79MeanProperties & prop,const int & cycle){
        db.store<PRE79MeanProperties>(prop,boostbase::where
          (id_kw,id)
          (user_label,user)
          (L_label,settings.lattice.L)
          (W_label,settings.lattice.W)
          (H_label,settings.lattice.H)
          (cycle_label,cycle)
          (ncycles_label,ncycles)
          (temperature_label,settings.hamiltonian.temperature)
          (lambda_label,settings.hamiltonian.lambda)
          (tau_label,settings.hamiltonian.tau)
          (type_label,current_properties_kw)
        );
    }
    void StoreFinalProperties(const Settings & settings,PRE79MeanProperties & prop){
        db.store<PRE79MeanProperties>(prop,boostbase::where
          (id_kw,id)
          (user_label,user)
          (L_label,settings.lattice.L)
          (W_label,settings.lattice.W)
          (H_label,settings.lattice.H)
          (cycle_label,ncycles-1)
          (ncycles_label,ncycles)
          (temperature_label,settings.hamiltonian.temperature)
          (lambda_label,settings.hamiltonian.lambda)
          (tau_label,settings.hamiltonian.tau)
          (type_label,final_properties_kw)
        );
    }
    boostbase::base & GetDB(){
        return db;
    }
    
};


#endif	/* _SIMULATIONDB_H */


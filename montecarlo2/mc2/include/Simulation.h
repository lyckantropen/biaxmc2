/* 
 * File:   Simulation.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 19:32
 */

#ifndef _SIMULATION_H
#define	_SIMULATION_H

#include "Metropolis.h"
#include "Lattice.h"
#include "std.h"

class Simulation {
    int acc_idx;
    int ncycles;
protected:
    virtual void DoIterate() {}
    Simulation(const int & nc){
        ncycles=nc;
        acc_idx=-1;
    }
public:
    bool Iterate(){
        if((acc_idx+1)>=ncycles) return false;
        acc_idx++;
        DoIterate();
        return true;
    }

    //Accessors
    const int & GetNCycles() const {
        return ncycles;
    }
    const int & GetAccIdx() const {
        return acc_idx;
    }
};

class LatticeSimulation:public Simulation {
    Lattice * lat;
    Hamiltonian * H;
    Metropolis * metropolis;
    virtual void DoIterate(){
        lat->Sweep(metropolis);
    }
public:
    LatticeSimulation(Hamiltonian * h=NULL,Lattice * l=NULL,Metropolis * metro=NULL,int nc=0):
    Simulation(nc),
    metropolis(metro),
    H(h),
    lat(l)
    {}
};

#endif	/* _SIMULATION_H */


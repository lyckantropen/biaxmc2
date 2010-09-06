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
#include "valarray_external.h"
#include "Statistical.h"

class Simulation {
    int acc_idx;
    int ncycles;
protected:
    virtual void DoIterate() {}
    Simulation(const int & nc,const int & sc=0){
        ncycles=nc;
        acc_idx=sc-1;
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
    double  accepted_fraction;
    vect acf;
    virtual void DoIterate(){
        accepted_fraction = double(lat->Sweep(metropolis))/double(lat->GetN());
        if(accepted_fraction>metropolis->GetAccULimit() || accepted_fraction<metropolis->GetAccLLimit())
            metropolis->AdjustRadius(lat);
        acf[GetAccIdx()]=accepted_fraction;
    }
public:
    LatticeSimulation(Hamiltonian * h=NULL,Lattice * l=NULL,Metropolis * metro=NULL,int nc=0,int startc=0):
    Simulation(nc,startc),
    metropolis(metro),
    H(h),
    lat(l),
    accepted_fraction(1.0)
    {
        acf.resize(nc,0.0);
    }

    const double & GetAcceptance() const {
        return accepted_fraction;
    }

    const double & GetMeanAcceptance() const {
        return Mean(acf,0,GetAccIdx()+1);
    }
};

#endif	/* _SIMULATION_H */


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
    Simulation(const int & nc=0,const int & sc=0){
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
    shared_ptr<Lattice> lat;
    shared_ptr<Hamiltonian>  H;
    shared_ptr<Metropolis> metropolis;
    double  accepted_fraction;
    double  accepted_fraction_p;
    vect acf;
    vect acfp;
    virtual void DoIterate(){
        int acc_rot=0,acc_p=0;
        lat->Sweep(metropolis,acc_rot,acc_p);
        accepted_fraction = double(acc_rot)/double(lat->GetN());
        accepted_fraction_p = double(acc_p)/double(lat->GetN());

        if(accepted_fraction>metropolis->GetAccULimit() || accepted_fraction<metropolis->GetAccLLimit())
            metropolis->AdjustRadius(lat);

        //to nie ma sensu!
        //if(accepted_fraction_p>metropolis->GetAccULimit() || accepted_fraction_p<metropolis->GetAccLLimit())
        //    metropolis->AdjustParityProb(lat);

        acf[GetAccIdx()]=accepted_fraction;
        acfp[GetAccIdx()]=accepted_fraction_p;
    }
public:
    LatticeSimulation(shared_ptr<Hamiltonian> h=shared_ptr<Hamiltonian>(),shared_ptr<Lattice> l=shared_ptr<Lattice>(),shared_ptr<Metropolis> metro=shared_ptr<Metropolis>(),int nc=0,int startc=0):
    Simulation(nc,startc),
    metropolis(metro),
    H(h),
    lat(l),
    accepted_fraction(1.0)
    {
        acf.resize(nc,0.0);
        acfp.resize(nc,0.0);
    }/*
    LatticeSimulation(){
        lat=NULL;
        H=NULL;
        metropolis=NULL;
    }*/
    
    LatticeSimulation(const LatticeSimulation & s){
        lat=s.lat;
        H=s.H;
        metropolis=s.metropolis;
        accepted_fraction=s.accepted_fraction;
        accepted_fraction_p=s.accepted_fraction_p;
        acf.resize(s.acf.size(),0.0);
        acfp.resize(s.acfp.size(),0.0);
        acf=s.acf;
        acfp=s.acfp;
    }
    const LatticeSimulation & operator=(const LatticeSimulation & s){
        lat=s.lat;
        H=s.H;
        metropolis=s.metropolis;
        accepted_fraction=s.accepted_fraction;
        accepted_fraction_p=s.accepted_fraction_p;
        acf.resize(s.acf.size(),0.0);
        acfp.resize(s.acfp.size(),0.0);
        acf=s.acf;
        acfp=s.acfp;
        return *this;
    }

    const double & GetAcceptance() const {
        return accepted_fraction;
    }
    const double & GetAcceptanceP() const {
        return accepted_fraction_p;
    }

    const double & GetMeanAcceptance() const {
        return Mean(acf,0,GetAccIdx()+1);
    }
    const double & GetMeanAcceptanceP() const {
        return Mean(acfp,0,GetAccIdx()+1);
    }
    void SetLattice(shared_ptr<Lattice> l){
        lat=l;
    }
};

#endif	/* _SIMULATION_H */


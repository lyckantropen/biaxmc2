/*
 * File:   Simulation.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 19:32
 */

#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "Metropolis.h"
#include "Lattice.h"
#include "std.h"
#include "Maths.h"
#include "Statistical.h"

///
/// A Monte Carlo simulation prototype
///

class Simulation
{
    int acc_idx;
    int ncycles;
protected:
    virtual void DoIterate() {}
    Simulation(const int & nc = 0, const int & sc = 0)
    {
        ncycles = nc;
        acc_idx = sc - 1;
    }
public:
    bool Iterate()
    {
        if((acc_idx + 1) >= ncycles) return false;
        acc_idx++;
        DoIterate();
        return true;
    }

    //Accessors
    const int & GetNCycles() const
    {
        return ncycles;
    }
    const int & GetAccIdx() const
    {
        return acc_idx;
    }
};

///
/// A general lattice simulation class. In principle, you
/// could isolate this class and build your own program
/// by implementing a different Hamiltonian and Metropolis
/// algorightm.
///
class LatticeSimulation: public Simulation
{
    Lattice * lat;
    Hamiltonian *  H;
    Metropolis * metropolis;
    double  accepted_fraction;
    double  accepted_fraction_p;
    vect acf;
    vect acfp;
    virtual void DoIterate()
    {
        int acc_rot = 0, acc_p = 0;
        lat->Sweep(metropolis, acc_rot, acc_p);
        accepted_fraction = double(acc_rot) / double(lat->GetN());
        accepted_fraction_p = double(acc_p) / double(lat->GetN());

        if(accepted_fraction > metropolis->GetAccULimit() || accepted_fraction < metropolis->GetAccLLimit())
            metropolis->AdjustRadius(lat);

        //to nie ma sensu!
        //if(accepted_fraction_p>metropolis->GetAccULimit() || accepted_fraction_p<metropolis->GetAccLLimit())
        //    metropolis->AdjustParityProb(lat);

        acf[GetAccIdx()] = accepted_fraction;
        acfp[GetAccIdx()] = accepted_fraction_p;
    }
public:
    LatticeSimulation(Hamiltonian * h = NULL, Lattice * l = NULL, Metropolis * metro = NULL, int nc = 0, int startc = 0):
        Simulation(nc, startc),
        metropolis(metro),
        H(h),
        lat(l),
        accepted_fraction(1.0)
    {
        acf.resize(nc, 0.0);
        acfp.resize(nc, 0.0);
    }

    const double & GetAcceptance() const
    {
        return accepted_fraction;
    }
    const double & GetAcceptanceP() const
    {
        return accepted_fraction_p;
    }

    const double & GetMeanAcceptance() const
    {
        return Mean(acf, 0, GetAccIdx() + 1);
    }
    const double & GetMeanAcceptanceP() const
    {
        return Mean(acfp, 0, GetAccIdx() + 1);
    }
    void SetLattice(Lattice * l)
    {
        lat = l;
    }
};

#endif  /* _SIMULATION_H */


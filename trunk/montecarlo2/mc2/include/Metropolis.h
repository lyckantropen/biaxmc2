/*
 * File:   Metropolis.h
 * Author: karol
 *
 * Created on 18 listopad 2009, 00:12
 */

#ifndef _METROPOLIS_H
#define _METROPOLIS_H

#include "MCProto.h"
#include "ILoggable.h"
#include "Settings.h"

/**
 * Implementacja algorytmu Metropolisa dla naszego przypadku
 * Cykl Monte Carlo wykonywany jest z użyciem (przykładowo) tej klasy
 * bezpośrednio przez klasę Lattice i, następnie pośrednio, Particle.
 */
class Hamiltonian;
class Lattice;
class Metropolis: public MCProto, protected ILoggable
{
    const Hamiltonian         * hamiltonian;
    double              radius; ///<promień dający odpowiednią akceptację ruchów
    double  acc_llimit, acc_ulimit;
    double  parity_prob;
    const Settings * settings;
public:
    Metropolis(const Settings & set, const Hamiltonian * h);
    Metropolis(const Metropolis & s);
    const Metropolis & operator=(const Metropolis & s);
    virtual vect OrientationNudge(const vect & old) const;
    virtual short ParityNudge(const short & old) const;
    virtual bool Accept(const double & dE)const;
    virtual const Hamiltonian * GetHamiltonian() const;
    double MeasureAccepted(Lattice * lat);
    void AdjustRadius(Lattice * lat, const double & decimation = 0.01);

    virtual void SetStream(std::ostream * os);
    const double & GetAccLLimit() const;
    const double & GetAccULimit() const;
    const double & GetParityProb() const;
};


#endif  /* _METROPOLIS_H */


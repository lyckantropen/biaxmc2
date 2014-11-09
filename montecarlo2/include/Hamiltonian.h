/*
 * File:   hamiltonian.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 16:11
 */

#ifndef _HAMILTONIAN_H
#define _HAMILTONIAN_H

class Particle;

/// purely virtual class, base class for any two-particle hamiltonian
class Hamiltonian
{
public:
    virtual double TwoParticleEnergy(const Particle &, const Particle &) const = 0;
    virtual double ExternalInteractionEnergy(const Particle &) const = 0;

};

#endif  /* _HAMILTONIAN_H */


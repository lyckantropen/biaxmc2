/* 
 * File:   hamiltonian.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 16:11
 */

#ifndef _HAMILTONIAN_H
#define	_HAMILTONIAN_H

class Particle;

class Hamiltonian {
public:
    virtual double TwoParticleEnergy(const Particle &, const Particle &)=0;
    virtual double ExternalInteractionEnergy(const Particle &)=0;

};

#endif	/* _HAMILTONIAN_H */


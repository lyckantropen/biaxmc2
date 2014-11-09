/*
 * File:   mcproto.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 16:50
 */

#ifndef _MCPROTO_H
#define _MCPROTO_H

#include "Maths.h"
#include "Hamiltonian.h"

/// Prototype for a Monte Carlo mechanism, such as Metropolis
class MCProto
{
public:
    ///najbardziej ogólnie, krok nastepny zależy od poprzedniego stanu
    virtual vect OrientationNudge(const vect & old) const = 0;
    ///ogólnie może zależeć od starej parzystości, ale zwykle tak nie jest
    virtual short ParityNudge(const short & old) const = 0;
    ///ogół decyzji prowadzących do zaakceptowania ruchu, na podstawie różnicy w energii
    virtual bool Accept(const double & dE)const = 0;
    virtual const Hamiltonian * GetHamiltonian()const = 0;
};


#endif  /* _MCPROTO_H */


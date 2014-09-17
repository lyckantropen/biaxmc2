/*
 * File:   AutoCorrelationTimeCalculator.h
 * Author: karolnew
 *
 * Created on 30 listopad 2011, 11:22
 */

#ifndef AUTOCORRELATIONTIMECALCULATOR_H
#define AUTOCORRELATIONTIMECALCULATOR_H

#include "PRE79StandardHamiltonian.h"
#include "Lattice.h"
#include "Statistical.h"
#include "Contractions.h"
#include "PRE79SpatialCorrelations.h"

#include "serializer.h"
#include "ILoggable.h"

class AutoCorrelationTimeCalculator: public ILoggable
{
    const Lattice * lattice;
    int ncycles;
    int acc_idx;
    vect ec;
    vect R;
    int t;


    Value CalculateMeanEnergy();
public:
    AutoCorrelationTimeCalculator(const Lattice * l, int _nc, int _nh);

    AutoCorrelationTimeCalculator(const AutoCorrelationTimeCalculator & s);

    const AutoCorrelationTimeCalculator & operator=(const AutoCorrelationTimeCalculator & s);
    void Update();
    const int & GetT() const;
    const vect & GetR() const;

};



#endif  /* AUTOCORRELATIONTIMECALCULATOR_H */


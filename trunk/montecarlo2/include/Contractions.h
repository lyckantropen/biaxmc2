/* 
 * File:   Contractions.h
 * Author: karol
 *
 * Created on 24 listopad 2009, 19:37
 */

#ifndef _CONTRACTIONS_H
#define	_CONTRACTIONS_H

#include "Maths.h"

/// Contraction over all indices of two 3-dimensional rank 2 symmetric tensors
inline double  Rank2Contraction(const vect & a, const vect & b)
{
    vect coeff(6);
    coeff[0]=coeff[3]=coeff[5]=1.0;
    coeff[1]=coeff[2]=coeff[4]=2.0;
    coeff*=a*b;
    return coeff.sum();
}
/// Contraction over all indices of two 3-dimensional rank 3 symmetric tensors
inline double  Rank3Contraction(const vect & a, const vect & b)
{
    vect coeff(10);
    coeff[0]=1;
    coeff[1]=3;
    coeff[2]=3;
    coeff[3]=1;
    coeff[4]=3;
    coeff[5]=6;
    coeff[6]=3;
    coeff[7]=3;
    coeff[8]=3;
    coeff[9]=1;
    coeff*=a*b;
    return coeff.sum();
}

#endif	/* _CONTRACTIONS_H */


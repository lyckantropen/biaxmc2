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
double  Rank2Contraction(const vect & a, const vect & b);
/// Contraction over all indices of two 3-dimensional rank 3 symmetric tensors
double  Rank3Contraction(const vect & a, const vect & b);

#endif	/* _CONTRACTIONS_H */


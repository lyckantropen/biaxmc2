/* 
 * File:   Contractions.h
 * Author: karol
 *
 * Created on 24 listopad 2009, 19:37
 */

#ifndef _CONTRACTIONS_H
#define	_CONTRACTIONS_H

#include "valarray_external.h"

// kontrakcja tensora 3x3 kodowanego 6-cioma składowymi
extern double  Rank2Contraction(const vect & a, const vect & b);
// kantrakcja tensora 3x3x kodowanego 10-cioma składowymi
extern double  Rank3Contraction(const vect & a, const vect & b);

#endif	/* _CONTRACTIONS_H */


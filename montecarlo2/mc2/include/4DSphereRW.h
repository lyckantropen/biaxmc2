/* 
 * File:   4DSphereRW.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 13:16
 */

#ifndef _4DSPHERERW_H
#define	_4DSPHERERW_H

#include "Random.h"
#include "Maths.h"

/// Returns a random point on a 4-dimensional sphere of
/// radius 'r'
/// (this algorithm is due to Marsaglia)
vect RandomPointOn4DSphereMarsaglia(const double & r);

/// Returns a random point on a 4-dimensional sphere of
/// radius 'r'
/// (an older and less efficient version)
vect RandomPointOn4DSphereOld(const double & r);

/// Given a point on a unit 4-sphere 'x', returns a random
/// point on a unit 4-sphere at distance 'r' from 'x'.
///
/// This is the basic move of a random walk on a 4-dimensional
/// sphere, which is used to generate random incremental
/// rotations by explointing the quaternion algebra.
vect RandomWalkOn4DSphere(const double & r,const vect & x);


#endif	/* _4DSPHERERW_H */


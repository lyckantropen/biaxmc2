/*
 * File:   Maths.h
 * Author: karol
 *
 * Created on 16 listopad 2009, 17:55
 *
 * Various helper mathematical functions and typedefs
 *
 */

#ifndef _Maths_H
#define _Maths_H

#include <valarray>
#include <iostream>
#include <vector>
#include "serializer.h"

///
/// Throughout biaxmc2 we use the excellent valarray template from
/// STL for computations on arrays. It defines numerous operators
/// and functions which can be performed over the entire array
/// and is optimized for speed.
///
/// Note on matrices: All symmetric 3x3 matrices are expressed
/// as 6-component 'vect' arrays as such:
///     / v[0] v[1] v[4] \
/// X = | v[1] v[2] v[3] |
///     \ v[4] v[5] v[6] /
///
typedef std::valarray<double>   vect;

/// Just the sign
#define sgn(x) (std::ceil((std::abs(x))/(x)))

/** Several handy functions for vect **/

inline double DotProduct(const vect & a, const vect & b)
{
    return (a*b).sum();
}

///euclidean norm
inline double Norm(const vect & v)
{
    vect result = std::pow(v,2.0);
    return std::sqrt(result.sum());
}

/// Identity matrix
/// Note: only dimension 3 is supported
inline vect Identity(const int & dim = 3)
{
    //TODO: wyższe wymiary niż 3
    vect out(6);
    out[0]=out[3]=out[5]=1.0;
    out[1]=out[2]=out[4]=0.0;
    return out;
}

/// Matrix dot product for 2 SYMMETRIC matrices parametrized
/// as 6-component 'vect' arrays
inline double MatrixDotProduct(const vect & a, const vect & b)
{
    return a[0]*b[0]+
            a[3]*b[3]+
            a[5]*b[5]+
            2.0*a[1]*b[1]+
            2.0*a[2]*b[2]+
            2.0*a[4]*b[4];
}


extern double Minimum(const vect & a);
extern int MinimumIndex(const vect & a);

/// Prints the array in a Mathematica form: { vect[0], vect[1], ... }
extern std::string MathematicaForm(const vect &);
extern std::ostream & operator<<(std::ostream & o, const vect & v);

/// This is a masking operator for the 'vect' type needed for some operations
template<class V>
std::vector<V>  operator|(const std::vector<V> & src, const std::valarray<bool> & mask)
{
    std::vector<V> ret;
    for(int i = 0; i < src.size(); i++)
        if(mask[i]) ret.push_back(src[i]);
    return ret;
}

#endif  /* _Maths_H */


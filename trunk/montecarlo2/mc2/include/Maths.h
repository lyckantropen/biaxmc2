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

extern double  DotProduct(const vect & a, const vect & b);

///euclidean norm
extern double  Norm(const vect &);

/// Identity matrix
/// Note: only dimension 3 is supported
extern vect Identity(const int & dim = 3);

/// Matrix dot product for 2 SYMMETRIC matrices parametrized
/// as 6-component 'vect' arrays
extern double MatrixDotProduct(const vect & a, const vect & b);


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

/// Reimplementation of serialization operators needed for GCC
//#ifdef
//template<class stream_t>
//void operator|(boostbase::outserializer<stream_t> & s, std::vector<vect> & v)
//{
//    s.template operator | <vect>(v);
//}
//template<class stream_t>
//void operator|(boostbase::inserializer<stream_t> & s, std::vector<vect> & v)
//{
//    s.template operator | <vect>(v);
//}
//#endif

#endif  /* _Maths_H */


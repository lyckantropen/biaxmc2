/* 
 * File:   valarray_external.h
 * Author: karol
 *
 * Created on 16 listopad 2009, 17:55
 */

#ifndef _VALARRAY_EXTERNAL_H
#define	_VALARRAY_EXTERNAL_H

#include <valarray>
#include <iostream>
#include "serializer.h"
typedef std::valarray<double>   vect;

extern double  DotProduct(const vect & a, const vect & b);
extern std::ostream & operator<<(std::ostream & o, const vect & v);
//norma euklidesowa
extern double  Norm(const vect &);
//macierz jednostkowa
extern vect Identity(const int & dim);
extern std::string MathematicaForm(const vect &);
extern double MatrixDotProduct(const vect & a, const vect & b);
extern double Minimum(const vect & a);
extern int MinimumIndex(const vect & a);

#define sgn(x) (std::ceil((std::abs(x))/(x)))

#ifndef __ICC
template<class stream_t>
void operator|(boostbase::outserializer<stream_t> & s, std::vector<vect> & v){
    s.operator|<vect>(v);
}
template<class stream_t>
void operator|(boostbase::inserializer<stream_t> & s, std::vector<vect> & v){
    s.operator|<vect>(v);
}
#endif

#endif	/* _VALARRAY_EXTERNAL_H */


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

template<class stream_t>
void operator|(boostbase::outserializer<stream_t> & s, std::vector<vect> & v){
    s.operator|<vect>(v);
}
template<class stream_t>
void operator|(boostbase::inserializer<stream_t> & s, std::vector<vect> & v){
    s.operator|<vect>(v);
}

#endif	/* _VALARRAY_EXTERNAL_H */


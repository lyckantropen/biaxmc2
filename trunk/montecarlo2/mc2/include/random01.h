/* 
 * File:   random01.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 12:50
 */

#ifndef _RANDOM01_H
#define	_RANDOM01_H
#include "boost.h"
#include "std.h"
#include "singleton.h"
#include "valarray_external.h"

class rangen:public Singleton<rangen> {
    boost::mt19937  rng;
    boost::uniform_01<boost::mt19937>   uni01;
    //boost::minstd_rand   rng;
    //boost::uniform_01<boost::minstd_rand>    uni01;
    friend class Singleton<rangen>;
protected:
    rangen():uni01(rng){}
public:
    double Gen01(){
        return uni01();
    }
};

extern double random01();
extern int plusminusone();
extern vect    RandomPointOn4DSphereMarsaglia(const double & r);
extern vect    RandomPointOn4DSphereOld(const double & r);


#endif	/* _RANDOM01_H */


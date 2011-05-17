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
#include "WydroRNG.h"

/*
 * Generator random01 jest głównym generatorem liczb losowych, ale zastrzegamy, żeby istniała osobna instancja
 * dla każdego wątku. Każda instancja jest inicjalizowana z generatora, który jest wspólny dla wszystkich wątków.
 */

//class rangen;
//extern rangen * rg;
//#pragma omp threadprivate(rg)
///ten generator jest wspólny dla każdego wątku i służy tylko do inicjalizacji instancji generatora rg (dyrektywa shared)
extern boost::rand48   rng2;


class rangen {
    boost::mt19937  rng;
    boost::uniform_01<boost::mt19937>   uni01;
    //boost::minstd_rand   rng;
    //boost::uniform_01<boost::minstd_rand>    uni01;
    //friend class Singleton<rangen>;
public:
    rangen():rng(rng2()),uni01(rng){}
    rangen(const rangen & r):rng(rng2()),uni01(rng){
    }
    double Gen01(){
        return uni01();
    }
    double operator()(){
        return uni01();
    }
};

class rangen_wd {
public:
    rangen_wd(){
        initializeRNG(seed1);
    }
    double Gen01(){
        return drandom(seed2);
    }
    double operator()(){
        return Gen01();
    }
};

///ten generator będzie osobny dla każdego wątku (dyrektywa threadprivate)
//extern rangen_wd random01;
extern rangen random01;

//extern double random01();
extern int plusminusone();
extern vect    RandomPointOn4DSphereMarsaglia(const double & r);
extern vect    RandomPointOn4DSphereOld(const double & r);


#endif	/* _RANDOM01_H */


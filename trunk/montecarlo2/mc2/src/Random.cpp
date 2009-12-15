/* 
 * File:   random01.cpp
 * Author: karol
 * 
 * Created on 17 listopad 2009, 12:50
 */

#include "random01.h"

double random01(){
    double retval = 0.0;
    // the rng is a singleton, so we need to make sure two threads don't access it simultaneously
    #pragma omp critical
    retval = rangen::Instance()->Gen01();
    return retval;
}
int plusminusone(){
    if(random01()>0.5)
        return 1;
    else return -1;
}
///this algorithm is due to Marsaglia
vect    RandomPointOn4DSphereMarsaglia(const double & r){
    vect result(4);
    double r1=0,r2=0;
    double y1=0,y2=0,y3=0,y4=0;
    do {
        y1=(1.0-2.0*random01());
        y2=(1.0-2.0*random01());
        r1=y1*y1+y2*y2;
    } while(r1>1);
    do {
        y3=(1.0-2.0*random01());
        y4=(1.0-2.0*random01());
        r2=y3*y3+y4*y4;
    } while(r2>1);
    double sr=std::sqrt((1-r1)/r2);

    result[0]=r*y1;
    result[1]=r*y2;
    result[2]=r*y3*sr;
    result[3]=r*y4*sr;
    return result;
}

vect    RandomPointOn4DSphereOld(const double & r){
	vect x(4);
	double r1;
	do
	{
		x[0]=(1-2*random01());
		x[1]=(1-2*random01());
		x[2]=(1-2*random01());
		x[3]=(1-2*random01());
		r1=Norm(x);
	} while (r1>1.0);
	x[0]=x[0]/r1*r;
	x[1]=x[1]/r1*r;
	x[2]=x[2]/r1*r;
	x[3]=x[3]/r1*r;
	return x;
}

template<> rangen * Singleton<rangen>::instance=NULL;

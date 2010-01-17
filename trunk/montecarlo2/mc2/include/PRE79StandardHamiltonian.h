/* 
 * File:   StandardHamiltonian.h
 * Author: karol
 *
 * Created on 18 listopad 2009, 01:03
 */

#ifndef _STANDARDHAMILTONIAN_H
#define	_STANDARDHAMILTONIAN_H

#include "Hamiltonian.h"
#include "Contractions.h"

class PRE79StandardHamiltonian:public Hamiltonian {
    double temperature ;
    double lambda;
    double tau;
    double h;   ///<second order field coupling
    double v1,v2,v3,vt,epsilon;
    double hxx,hyy,hzz;
public:
    PRE79StandardHamiltonian(const double & temp=1.0, const double & lam=0.0, const double & t=0.0, const double & hval=0.0){
        temperature=temp;
        tau=t;
        lambda=lam;
        h=hval;

        epsilon = 1.0/temperature;

        v1=epsilon*(2*lambda*lambda+sqrt6*lambda);
        v2=epsilon*(2*lambda*lambda-sqrt6*lambda);
        v3=epsilon*(3.0/2-lambda*lambda);
        vt=epsilon*tau;
        hxx = (lambda-1/sqrt6);
        hyy = -(lambda+1/sqrt6);
        hzz = std::sqrt(2/3);

    }
    virtual double ExternalInteractionEnergy(const Particle & p){
        // pole skierowane wzdłuż osi x
        return - (h/std::abs(h)) *h*h*(hxx*p.GetQX()[0]+hyy*p.GetQY()[0]+hzz*p.GetQZ()[0]);
        //return 0.0;
    }
    virtual double TwoParticleEnergy(const Particle & a, const Particle & b){
        double xx = DotProduct(a.GetEX(),b.GetEX());
        double yy = DotProduct(a.GetEY(),b.GetEY());
        double zz = DotProduct(a.GetEZ(),b.GetEZ());

        // uwaga - TEMPERATURA siedzi w stałych
        return - v1*xx*xx - v2*yy*yy - v3*zz*zz - vt*Rank3Contraction(a.GetT(),b.GetT());
 
    }

    //Accessors
    const double & GetTemperature() const {
        return temperature;
    }
    const double & GetLambda() const {
        return lambda;
    }
    const double & GetTau() const {
        return tau;
    }
};


#endif	/* _STANDARDHAMILTONIAN_H */


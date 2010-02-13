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

        //wartości Kacpra. Nie wiem skąd je wziął.
        //hxx = (lambda-1/sqrt6);
        //hyy = -(lambda+1/sqrt6);
        //hzz = std::sqrt(2/3);
        hxx = lambda;
        hyy = -lambda;
        hzz = std::sqrt(2./3.);

    }
    virtual double ExternalInteractionEnergy(const Particle & p){
        // pole skierowane wzdłuż osi z
        // uwaga - h jest równoważne h kwadrat w literaturze, a h/abs(h) jest równoważne delta epsilon i jest równe sgn(h)
        // tensor Q zawiera jeszcze część sferycznie symetryczną, którą pomijamy, choć można ją dodać (wynosi ona -1/sqrt(6))
        if(h==0.0) return 0.0;
        return - epsilon*(h/std::abs(h))*h*(hxx*p.GetQX()[5] + hyy*p.GetQY()[5] + hzz*p.GetQZ()[5] - 1/sqrt6);

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
    const double & GetH() const {
        return h;
    }
};


#endif	/* _STANDARDHAMILTONIAN_H */


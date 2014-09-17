/*
 * File:   StandardHamiltonian.h
 * Author: karol
 *
 * Created on 18 listopad 2009, 01:03
 */

#ifndef _STANDARDHAMILTONIAN_H
#define _STANDARDHAMILTONIAN_H

#include "Hamiltonian.h"
#include "Contractions.h"
#include "Particle.h"

class StandardL2Hamiltonian: public Hamiltonian
{
    double temperature ;
    double lambda;
    double tau;
    double kappa;
    double h;   ///<second order field coupling
    double v1, v2, v3, vt, vk, epsilon;
    double hxx, hyy, hzz;

    double LatticeCoupling(const Particle & a, const Particle & b) const
    {
        double p1 = a.GetParity();
        double p2 = b.GetParity();
        if(p1!=p2) return 0.;
        vect r = b.GetR() - a.GetR();
        ///account for boundary effects
        if(Norm(r) > 1.01)
        {
            //std::cout << r << " changed to ";
            r = -r;
            r /= Norm(r);
            //std::cout << r << std::endl;
        }
        for(int k = 0; k < 3; k++)
        {
            if(r[k] < 0.1) r[k] = std::floor(r[k]);
            if(r[k] > 0.5) r[k] = std::ceil(r[k]);
        }
        //r[0]=0.0;
        //r[2]=0.0;
        vect q=std::sqrt(1.5)*(-Identity(3)/3.+a.GetQZ())+lambda*(a.GetQX()-a.GetQY());
        vect s=std::sqrt(1.5)*(-Identity(3)/3.+b.GetQZ())+lambda*(b.GetQX()-b.GetQY());
     
        return (p1+p2)*(q[4]*r[1]*s[1]+q[0]*r[2]*s[1]-q[3]*r[2]*s[1]-q[0]*r[1]*s[2]+q[5]*r[1]*s[2]-q[4]*r[2]*s[2]-q[4]*r[0]*s[3]+q[3]*r[0]*s[4]-q[5]*r[0]*s[4]+q[1]*(r[0]*s[2]+r[2]*(-s[0]+s[3])-r[1]*s[4])+q[2]*(-r[0]*s[1]+r[2]*s[4]+r[1]*(s[0]-s[5]))+q[4]*r[0]*s[5]);
    }


public:
    StandardL2Hamiltonian(const double & temp = 1.0, const double & lam = 0.0, const double & t = 0.0, const double & hval = 0.0, const double & kap = 0.0)
    {
        temperature = temp;
        tau = t;
        lambda = lam;
        h = hval;
        kappa = kap;

        epsilon = 1.0 / temperature;

        v1 = epsilon * (2 * lambda * lambda + sqrt6 * lambda);
        v2 = epsilon * (2 * lambda * lambda - sqrt6 * lambda);
        v3 = epsilon * (3.0 / 2 - lambda * lambda);
        vt = epsilon * tau;
        vk = epsilon * kappa;

        //wartości Kacpra. Te są dobre, trzeba rozpisać jedynkę
        hxx = (lambda - 1 / sqrt6);
        hyy = -(lambda + 1 / sqrt6);
        hzz = std::sqrt(2 / 3);
        //hxx = lambda;
        //hyy = -lambda;
        //hzz = std::sqrt(2./3.);

    }
    virtual double ExternalInteractionEnergy(const Particle & p) const
    {
        // pole skierowane wzdłuż osi z, gdy pole dodatnie, a wzdłuż osi x, gdy ujemne. Ma to uzasadnienie ze względu na stan podstawowy
        // uwaga - h jest równoważne h kwadrat w literaturze, a h/abs(h) jest równoważne delta epsilon i jest równe sgn(h)
        //int sign = std::ceil(h)/std::abs(std::ceil(h));
        if(h == 0.0) return 0.0;
        //if(h>=0.0)
        return - epsilon * h * (hxx * p.GetQX()[5] + hyy * p.GetQY()[5] + hzz * p.GetQZ()[5]);
        //else
        //  return - epsilon*h*(hxx*p.GetQX()[0] + hyy*p.GetQY()[0] + hzz*p.GetQZ()[0]);

        //return 0.0;
    }
    virtual double TwoParticleEnergy(const Particle & a, const Particle & b) const
    {
        double xx = DotProduct(a.GetEX(), b.GetEX());
        double yy = DotProduct(a.GetEY(), b.GetEY());
        double zz = DotProduct(a.GetEZ(), b.GetEZ());

        // uwaga - TEMPERATURA siedzi w stałych
        return - v1 * xx * xx - v2 * yy * yy - v3 * zz * zz - vt * Rank3Contraction(a.GetT(), b.GetT()) - vk * LatticeCoupling(a, b);

    }

    void SetTemperature(const double & t)
    {
        temperature = t;
        epsilon = 1.0 / temperature;

        v1 = epsilon * (2 * lambda * lambda + sqrt6 * lambda);
        v2 = epsilon * (2 * lambda * lambda - sqrt6 * lambda);
        v3 = epsilon * (3.0 / 2 - lambda * lambda);
        vt = epsilon * tau;
    }

    //Accessors
    const double & GetTemperature() const
    {
        return temperature;
    }
    const double & GetLambda() const
    {
        return lambda;
    }
    const double & GetTau() const
    {
        return tau;
    }
    const double & GetH() const
    {
        return h;
    }
    const double & GetKappa() const
    {
        return kappa;
    }
};


#endif  /* _STANDARDHAMILTONIAN_H */


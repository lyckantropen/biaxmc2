/* 
 * File:   Metropolis.h
 * Author: karol
 *
 * Created on 18 listopad 2009, 00:12
 */

#ifndef _METROPOLIS_H
#define	_METROPOLIS_H

#include "MCProto.h"
#include "Hamiltonian.h"
#include "4DSphereRW.h"
#include "Lattice.h"
#include "ILoggable.h"
#include "Settings.h"

/**
 * Implementacja algorytmu Metropolisa dla naszego przypadku
 * Cykl Monte Carlo wykonywany jest z użyciem (przykładowo) tej klasy
 * bezpośrednio przez klasę Lattice i, następnie pośrednio, Particle.
 */
class Metropolis:public MCProto,protected ILoggable {
    Hamiltonian     *   hamiltonian;
    double              radius; ///<promień dający odpowiednią akceptację ruchów
    double  acc_llimit,acc_ulimit;
    const Settings & settings;
public:
    Metropolis(const Settings & set, Hamiltonian * h=NULL,const double & r=1):settings(set),hamiltonian(h),radius(r),
            acc_llimit(settings.simulation.metropolis_lower_acceptance_limit),
            acc_ulimit(settings.simulation.metropolis_higher_acceptance_limit){}
    virtual vect OrientationNudge(const vect & old){
        return RandomWalkOn4DSphere(radius,old);
    }
    virtual short ParityNudge(const short & old){
        
        if(random01()<settings.simulation.parity_flip_probability)
            return -old;
        else return old;
        
        //return plusminusone();
    }
    virtual bool Accept(const double & dE){
        if(dE<0)
            return true;
        else {
            double acceptance = std::exp(-dE);
            if(random01()<acceptance)
                return true;
            else
                return false;
        }

    }
    virtual Hamiltonian * GetHamiltonian(){
        return hamiltonian;
    }
    double MeasureAccepted(Lattice * lat){
        if(lat==NULL) return -1 ;
        double N=0.0;
        int acc_moves=0;
        Lattice testlat=*lat;

        double tries=5;
        N=testlat.GetN();
        for(int i=0;i<tries;i++){
            acc_moves+=testlat.Sweep(this);
        }
        return double(acc_moves)/(N*tries);

    }
    void AdjustRadius(Lattice * lat, const double & decimation=0.02){
        if(lat==NULL) return ;
        double acc_fraction=0.0;
        double N=0.0;
        int acc_moves=0;


        while(acc_fraction<acc_llimit || acc_fraction>acc_ulimit){
            Lattice testlat=*lat;
            acc_moves=0;
            double tries=2;
            N=testlat.GetN();
            for(int i=0;i<tries;i++){
                acc_moves+=testlat.Sweep(this);
            }
            acc_fraction = double(acc_moves)/(N*tries);
            
            
            if(acc_fraction<acc_llimit){
                radius*=(1.0-decimation);
                //Log() << "decimating down to r=" << radius << " because acc_fraction=" << acc_fraction << std::endl;
            }
            if(acc_fraction>acc_ulimit){
                radius*=(1.0+decimation);
                //Log() << "decimating up to r=" << radius << " because acc_fraction=" << acc_fraction << std::endl;
            }
            
            if(radius>=1.0){
                radius=0.999;
                break;
            }
            if(radius<=0.001){
                radius=0.001;
                break;
            }
            
        }
        //Log() << acc_fraction << " (" << acc_moves << "/" << N << ") accepted, radius " << radius << std::endl;
    }
    virtual void SetStream(std::ostream * os){
        ILoggable::SetStream(os);
    }
    const double & GetAccLLimit() const {
        return acc_llimit;
    }
    const double & GetAccULimit() const {
        return acc_ulimit;
    }
};


#endif	/* _METROPOLIS_H */


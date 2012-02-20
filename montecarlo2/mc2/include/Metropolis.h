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
class Metropolis:public MCProto,protected ILoggable,public boost::enable_shared_from_this<Metropolis> {
    shared_ptr<Hamiltonian>     hamiltonian;
    double              radius; ///<promień dający odpowiednią akceptację ruchów
    double  acc_llimit,acc_ulimit;
    double  parity_prob;
    Settings settings;
public:
    Metropolis(const Settings & set, shared_ptr<Hamiltonian> h= shared_ptr<Hamiltonian>(),const double & r=1):settings(set),hamiltonian(h),radius(set.simulation.radius),
            acc_llimit(set.simulation.metropolis_lower_acceptance_limit),
            acc_ulimit(set.simulation.metropolis_higher_acceptance_limit),
            parity_prob(set.simulation.parity_flip_probability)
    {}
    Metropolis(const Metropolis & s){
        hamiltonian=s.hamiltonian;
        radius=s.radius;
        acc_llimit=s.acc_llimit;
        acc_ulimit=s.acc_ulimit;
        parity_prob=s.parity_prob;
        settings=s.settings;
    }
    const Metropolis & operator=(const Metropolis & s){
        hamiltonian=s.hamiltonian;
        radius=s.radius;
        acc_llimit=s.acc_llimit;
        acc_ulimit=s.acc_ulimit;
        parity_prob=s.parity_prob;
        settings=s.settings;
        return *this;
    }
    virtual vect OrientationNudge(const vect & old) const {
        return RandomWalkOn4DSphere(radius,old);
    }
    virtual short ParityNudge(const short & old) const {
        
        if(random01()<parity_prob)
            return -old;
        else return old;
        
        //return plusminusone();
    }
    virtual bool Accept(const double & dE)const {
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
    virtual const shared_ptr<Hamiltonian> GetHamiltonian() const {
        return hamiltonian;
    }
    double MeasureAccepted(shared_ptr<Lattice> lat){
        if(lat==NULL) return -1 ;
        double N=0.0;
        int acc_moves=0;
        Lattice testlat=*lat;

        double tries=5;
        N=testlat.GetN();
        for(int i=0;i<tries;i++){
            int acc_rot=0,acc_p=0;
            testlat.Sweep(dynamic_pointer_cast<MCProto>(shared_from_this()),acc_rot,acc_p);
            acc_moves+=acc_p;
        }
        return double(acc_moves)/(N*tries);

    }
    void AdjustRadius(shared_ptr<Lattice> lat, const double & decimation=0.01){
        if(!settings.simulation.adjust_radius) return;

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
                int acc_rot=0,acc_p=0;
                testlat.Sweep(dynamic_pointer_cast<MCProto>(shared_from_this()),acc_rot,acc_p);
                acc_moves+=acc_rot;
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
        Log() << acc_fraction << " (" << acc_moves << "/" << N << ") accepted, radius " << radius << std::endl;
    }
    /*
     * to nie ma sensu! zawsze zaakceptowanych jest parity_prob*N ruchów!
    void AdjustParityProb(Lattice * lat, const double & decimation=0.02){
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
                int acc_rot=0,acc_p=0;
                testlat.Sweep(this,acc_rot,acc_p);
                acc_moves+=acc_p;
            }
            acc_fraction = double(acc_moves)/(N*tries);


            if(acc_fraction<acc_llimit){
                parity_prob*=(1.0+decimation);
                Log() << "decimating down to r=" << parity_prob << " because acc_fraction=" << acc_fraction << std::endl;
            }
            if(acc_fraction>acc_ulimit){
                parity_prob*=(1.0-decimation);
                Log() << "decimating up to r=" << parity_prob << " because acc_fraction=" << acc_fraction << std::endl;
            }

            if(parity_prob>=1.0){
                parity_prob=0.999;
                break;
            }
            if(parity_prob<=0.001){
                parity_prob=0.001;
                break;
            }

        }
        //Log() << acc_fraction << " (" << acc_moves << "/" << N << ") accepted, probability " << parity_prob << std::endl;
    }
     */
    virtual void SetStream(std::ostream * os){
        ILoggable::SetStream(os);
    }
    const double & GetAccLLimit() const {
        return acc_llimit;
    }
    const double & GetAccULimit() const {
        return acc_ulimit;
    }
    const double & GetParityProb() const {
        return parity_prob;
    }
};


#endif	/* _METROPOLIS_H */


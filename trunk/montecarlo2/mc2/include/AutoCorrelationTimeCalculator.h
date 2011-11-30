/* 
 * File:   AutoCorrelationTimeCalculator.h
 * Author: karolnew
 *
 * Created on 30 listopad 2011, 11:22
 */

#ifndef AUTOCORRELATIONTIMECALCULATOR_H
#define	AUTOCORRELATIONTIMECALCULATOR_H

#include "PRE79StandardHamiltonian.h"
#include "Lattice.h"
#include "Statistical.h"
#include "Contractions.h"
#include "PRE79SpatialCorrelations.h"

#include "serializer.h"
#include "ILoggable.h"

class AutoCorrelationTimeCalculator:public ILoggable {
    Lattice                     *lattice;
    int ncycles;
    int acc_idx;
    vect ec;
    vect R;
    
    Value CalculateMeanEnergy(){
        vect e(0.0,lattice->GetN());
        for(int site=0;site<lattice->GetN();site++){
            e[site]=lattice->GetParticles()[site].GetEnergy();
        }
        return Mean(e,0,0,e.size());
    }
public:
    AutoCorrelationTimeCalculator(Lattice * l, int _nc, int _nh):lattice(l),ncycles(_nc){
        ec.resize(_nc,0.0);
        R.resize(_nh,0.0);
        acc_idx=0;
        SetFile("autocorrelation");
    }
    void Update() {
        ec[acc_idx]=double(CalculateMeanEnergy());
        
        if(acc_idx>=(ncycles-1)) {
            Log() << "Autocorrelation function:\n" ;

            const int N = (ncycles-R.size()+1);
            double em = Mean(ec,0,0,ec.size());
            double emv = std::pow(Mean(ec,0,0,ec.size()).Error(),2);
            for(int k=0;k<R.size();k++){
                double q1=0,q2=0,q3=0;
                for(int i=0;i<N;i++){
                    q1+=ec[i]/N;
                    q2+=ec[i+k]/N;
                    q3+=(ec[i]-em)*(ec[i+k]-em)/N;
                }
                R[k]=q3/emv;
                //std::cout << q1 << ", " << q2 << ", " << q3 << std::endl;
                Log() << R[k] << std::endl;
            }
            
     
            acc_idx=0;
            return;
        }
        acc_idx++;
    }
        
};



#endif	/* AUTOCORRELATIONTIMECALCULATOR_H */


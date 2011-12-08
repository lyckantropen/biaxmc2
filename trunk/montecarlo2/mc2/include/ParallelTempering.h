/* 
 * File:   ParallelTempering.h
 * Author: karolnew
 *
 * Created on 3 grudzień 2011, 18:13
 */

#ifndef PARALLELTEMPERING_H
#define	PARALLELTEMPERING_H

#include "ILoggable.h"
//#include "PRE79Simulation.h"
//#include "PRE79StandardProperties.h"
#include "valarray_external.h"
#include "Statistical.h"
#include "omp.h"

class ParallelTempering: public ILoggable {
    //std::vector<LatticeSimulation>        simulations;
    std::vector<Lattice>                lattices;
    SimulationDB        database;
    vect beta;
    vect E;
    vect index;
    vect acc;
    int n;
    Settings settings;
public:
    ParallelTempering(const Settings & s):settings(s),database(s) {
        n = settings.openmp.number_of_threads; //(settings.scanning.end-settings.scanning.start)/settings.scanning.delta + 1;
        beta.resize(n,0.0);
        E.resize(n,0.0);
        acc.resize(n-1,0.0);
        index.resize(n,0);
        //simulations.resize(n,PRE79Simulation(settings));
        //TODO: inne stany początkowe
        if(settings.initial.isotropic && !settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(Lattice(settings.lattice.L,settings.lattice.W,settings.lattice.H,
                        Lattice::Isotropic));
        else if(settings.initial.isotropic && settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(Lattice(settings.lattice.L,settings.lattice.W,settings.lattice.H,
                        Lattice::IsotropicRighthanded));
        else if(settings.initial.biaxial && !settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(Lattice(settings.lattice.L,settings.lattice.W,settings.lattice.H,
                        Lattice::Biaxial));
        else if(settings.initial.biaxial && settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(Lattice(settings.lattice.L,settings.lattice.W,settings.lattice.H,
                        Lattice::BiaxialRighthanded));
        
            
            
        double delta = (1./settings.scanning.start - 1./settings.scanning.end)/n;

        for(int i=0;i<n;i++){
            index[i]=i;
            if(settings.scanning.separate_values)
                beta[i]=1./settings.scanning.values[i];
            else
                beta[i]=1./settings.scanning.end+i*delta;
        }
        SetStream(&std::cout);
        Log() << "Betas: " << beta << std::endl;
    }
    void Run() {
        
        omp_set_dynamic(0);
        omp_set_num_threads(n);
        
        int swapfq = settings.scanning.parallel_tempering_swapfq;
        //int accfq = 10*swapfq;
        
        //#pragma omp parallel for schedule(static) private(random01) ordered
        #pragma omp parallel private(random01)
        //for(int i=0;i<n;i++)
        {
            int i = omp_get_thread_num();
            {
            int nc = settings.simulation.thermalization_cycles;
            
            //#pragma omp critical
            //Log() << "Thread: " << i << std::endl;
            PRE79StandardHamiltonian H(1./beta[i],settings.hamiltonian.lambda,settings.hamiltonian.tau,settings.hamiltonian.h);
            AutoCorrelationTimeCalculator ac(&lattices[i],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length);
            PRE79StandardProperties prop(&lattices[i],nc/settings.simulation.measure_frequency);
            PRE79StandardProperties ptprop(&lattices[i],nc/swapfq);
            Metropolis m(settings,&H);
            m.AdjustRadius(&lattices[i]);
            LatticeSimulation simulation(&H,&lattices[i],&m,nc);
            
            //ac.SetStream(&std::cout);
            for(int k=0;k<nc;k++){
                //H.SetTemperature(1./beta[index[i]]);
                H.SetTemperature(1./beta[i]);
                simulation.SetLattice(&lattices[index[i]]);
                prop.SetLattice(&lattices[index[i]]);
                
                simulation.Iterate();
                ac.Update();
                if(k%20==0){
                        prop.Update(k/settings.simulation.measure_frequency,&H,&ac);
                        //E[i]=prop.EnergyEvolution()[prop.GetAccIdx()];
                }
                if(k>0 && k%swapfq==0 && i<(n-1))
                {          
                   ptprop.Update(k/swapfq,&H,&ac);
                   E[i]=ptprop.EnergyEvolution()[ptprop.GetAccIdx()];

                   #pragma omp critical
                   Swap(i+1,i,k/swapfq);
                   
                   Log() << "New indices: " << index << std::endl;
                   Log() << "Simulation " << i << " has Beta " << /*beta[index[i]]*/ beta[i] << " but index " << index[i] << std::endl;
                   Log() << "No of accepted swaps: " << acc << std::endl;
                   
                }
              
                //#pragma omp barrier
                /*
                #pragma omp master
                if(k>2*swapfq && k%swapfq==0 && i==(n-1))
                {
                    //#pragma omp barrier
                
                    double aptm=0.0;
                                      
                    for(int i=1;i<n;i++)
                        aptm+=acc[i-1]*(beta[i]-beta[i-1]);
                    
                    double lambda = (beta[n-1]-beta[0])/aptm;
                    
                    Log() << "lambda = " << lambda << std::endl;
                    
                    vect oldbeta = beta;
                    for(int i=1;i<n;i++){
                        beta[i]=beta[i-1]+lambda*acc[i-1]*(oldbeta[i]-oldbeta[i-1]);
                        if(beta[i]==beta[i-1]){
                            beta[i-1]*=0.95;
                            beta[i]*=1.05;
                        }
                    }
                    Log() << "New betas: " << beta << std::endl;
                }*/
                 
               
                
            }
            //H.SetTemperature(1./beta[index[i]]);
            H.SetTemperature(1./beta[i]);
            //dotermalizowanie
            
            Log() << "Supplementary thermalization of " << settings.simulation.supplementary_thermalization_cycles << " cycles.\n";
            for(int k=0;k<settings.simulation.supplementary_thermalization_cycles;k++)
                simulation.Iterate();
            
            PRE79MeanProperties mprop(prop,H);
            Settings cur = settings;
            cur.hamiltonian.temperature = H.GetTemperature();//1./beta[index[i]];
            database.StoreFinalLattice(cur,lattices[i]);
            database.StoreThermalizationHistory(cur,prop);
            database.StoreProperties(cur,mprop,nc-1);
            }
            
                    ////////////////////////
                    
                    ////////////////////////
            
            if(settings.simulation.production_cycles>0)
            {
            int nc = settings.simulation.production_cycles;
            Log() << "Production with " << settings.simulation.production_cycles << " cycles\n";
            //#pragma omp critical
            //Log() << "Thread: " << i << std::endl;
            PRE79StandardHamiltonian H(1./beta[i],settings.hamiltonian.lambda,settings.hamiltonian.tau,settings.hamiltonian.h);
            AutoCorrelationTimeCalculator ac(&lattices[i],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length);
            PRE79StandardProperties prop(&lattices[i],nc/settings.simulation.measure_frequency);
            PRE79StandardProperties ptprop(&lattices[i],nc/swapfq);
            Metropolis m(settings,&H);
            m.AdjustRadius(&lattices[i]);
            LatticeSimulation simulation(&H,&lattices[i],&m,nc);
            
            //ac.SetStream(&std::cout);
            for(int k=0;k<nc;k++){
                //H.SetTemperature(1./beta[index[i]]);
                H.SetTemperature(1./beta[i]);
                simulation.SetLattice(&lattices[index[i]]);
                prop.SetLattice(&lattices[index[i]]);
                
                simulation.Iterate();
                ac.Update();
                if(k%20==0){
                        prop.Update(k/settings.simulation.measure_frequency,&H,&ac);
                        //E[i]=prop.EnergyEvolution()[prop.GetAccIdx()];
                }
                if(k>0 && k%swapfq==0 && i<(n-1))
                {          
                   ptprop.Update(k/swapfq,&H,&ac);
                   E[i]=ptprop.EnergyEvolution()[ptprop.GetAccIdx()];

                   #pragma omp critical
                   Swap(i+1,i,k/swapfq);
                   
                   Log() << "New indices: " << index << std::endl;
                   Log() << "Simulation " << i << " has Beta " << /*beta[index[i]]*/ beta[i] << " but index " << index[i] << std::endl;
                   Log() << "No of accepted swaps: " << acc << std::endl;
                   
                }
              
                //#pragma omp barrier
                /*
                #pragma omp master
                if(k>2*swapfq && k%swapfq==0 && i==(n-1))
                {
                    //#pragma omp barrier
                
                    double aptm=0.0;
                                      
                    for(int i=1;i<n;i++)
                        aptm+=acc[i-1]*(beta[i]-beta[i-1]);
                    
                    double lambda = (beta[n-1]-beta[0])/aptm;
                    
                    Log() << "lambda = " << lambda << std::endl;
                    
                    vect oldbeta = beta;
                    for(int i=1;i<n;i++){
                        beta[i]=beta[i-1]+lambda*acc[i-1]*(oldbeta[i]-oldbeta[i-1]);
                        if(beta[i]==beta[i-1]){
                            beta[i-1]*=0.95;
                            beta[i]*=1.05;
                        }
                    }
                    Log() << "New betas: " << beta << std::endl;
                }*/
                 
               
                
            }
            //H.SetTemperature(1./beta[index[i]]);
            H.SetTemperature(1./beta[i]);
            
            Log() << "Saving final properties for beta " << beta[i] << std::endl;
            PRE79MeanProperties mprop(prop,H);
            Settings cur = settings;
            cur.hamiltonian.temperature = H.GetTemperature();//1./beta[index[i]];
            database.StoreFinalLattice(cur,lattices[i]);
            database.StorePropertiesEvolution(cur,prop);
            database.StoreFinalProperties(cur,mprop);
            }
            
        }
    }
    void Swap(const int & u1, const int & u2,const int &k) {
        double dB = beta[u2]-beta[u1];
        double dE = E[u1] - E[u2];
        double frac = std::exp(-dB*dE);
        if(frac>=1.0) frac=1.0;
        
        double delta = frac - acc[u2];
        acc[u2]+=delta/(k+1);
        //acc[u2]=frac;
        
        //std::cout << "dB=" << dB << " dE=" << dE << ", exp(-dBdE)=" << frac << std::endl;
        //acc[u2]+=frac;
        //delete simulations[index[k]];
        if(frac>=1.0 or random01()<frac){
            int ii = index[u1];
            index[u1]=index[u2];
            index[u2]=ii;
            //acc[u2]+=1.;
        }/*
        else if(random01()<frac){
            int ii = index[u1];
            index[u1]=index[u2];
            index[u2]=ii;
            //acc[u2]+=1.;
        }*/      
    }
    /*
    void Run() {
        //int ncycles = 100000;
        int swapfq = 20;
        int swaps = settings.simulation.thermalization_cycles/swapfq;
        
        omp_set_dynamic(0);
        omp_set_num_threads(n);
        
        //#pragma omp parallel for
        for(int F=0;F<swaps;F++)
        #pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++)
        {
            
            / *foreach(int & u,index){
                std::cout << u << ",";
            }
            std::cout << std::endl;
            std::cout << "i: " << i << ", F: " <<F << std::endl;* /
            
            Settings current_settings = settings;
            current_settings.hamiltonian.temperature = 1./beta[i];
            current_settings.simulation.thermalization_cycles = settings.simulation.thermalization_cycles/swaps;
            current_settings.simulation.supplementary_thermalization_cycles = current_settings.simulation.thermalization_cycles;
            current_settings.output.save_thermalization_properties = false;
            current_settings.simulation.adjust_radius = false;
            current_settings.simulation.calculate_time = false;
            current_settings.output.report_progress = false;
            current_settings.simulation.measure_acceptance = false;
            if(F>0){
                PRE79Simulation s(current_settings,lattices[index[i]]);
                //s.SetStream(&std::cout);
                lattices[index[i]] = s.Thermalize();
                E[i]=s.GetThermalizationProperties().EnergyEvolution()[current_settings.simulation.thermalization_cycles/10-1];
            }
            else{
                PRE79Simulation s(current_settings);
                //s.SetStream(&std::cout);
                lattices[index[i]] = s.Thermalize();
                E[i]=s.GetThermalizationProperties().EnergyEvolution()[current_settings.simulation.thermalization_cycles/10-1];
            }
            if(F==swaps-1){
                SimulationDB db(settings);
                db.StoreFinalLattice(current_settings,lattices[index[i]]);
            }
            

            #pragma omp critical
            {
                if(omp_get_thread_num()==0){
                for(int k=0;k<(n-1);k++){
                    //int idx = settings.simulation.thermalization_cycles/swaps-1;//simulations[index[k+1]].GetThermalizationProperties().GetAccIdx() ;
                    double dB = beta[k+1]-beta[k];
                    double dE = E[k+1] - E[k];
                    //delete simulations[index[k]];
                    if(dE<0){
                        int ii = index[k+1];
                        index[k+1]=index[k];
                        index[k]=ii;
                        acc[k]+=1.;
                    }
                    else if(std::exp(-dB*dE)>random01()){
                        int ii = index[k+1];
                        index[k+1]=index[k];
                        index[k]=ii;
                        acc[k]+=1.;
                    }
                    //std::cout << "k: " << k << ", dB: " << dB << ", dE: " << dE << std::endl;
                    
                }
                / *
                foreach(int & u,index){
                    std::cout << u << ",";
                }
                std::cout << std::endl;
                 * /
                }
            }

            
            
        }
        
    }*/
    
};

#endif	/* PARALLELTEMPERING_H */


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
    std::vector<LatticeSimulation*>      simulations;
    std::vector<Lattice>                lattices;
    std::vector<PRE79StandardProperties*> prop;
    std::vector<PRE79StandardHamiltonian*> H;
    std::vector<Metropolis*>     m;
    std::vector<AutoCorrelationTimeCalculator*> ac;
    
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
        Log() << "Beta distance: " << delta << std::endl;

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
    void Run2() {
        //liczba cykli termalizacyjnych
        int nc = settings.simulation.thermalization_cycles;

        
        omp_set_dynamic(0);
        omp_set_num_threads(n);
        
        for(int i=0;i<n;i++)
        {
            
            H.push_back(new PRE79StandardHamiltonian(1./beta[i],settings.hamiltonian.lambda,settings.hamiltonian.tau,settings.hamiltonian.h));
            m.push_back(new Metropolis(settings,H.back()));
            prop.push_back(new PRE79StandardProperties(&lattices[i],nc/settings.simulation.measure_frequency));
            simulations.push_back(new LatticeSimulation(H.back(),&lattices[i],m.back(),nc));
            ac.push_back(new AutoCorrelationTimeCalculator(&lattices[i],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length));
        }
        Log() << "Thermalization\n";
        InnerLoop(nc);
        
        Log() << "Saving thermalization properties\n";
        
        for(int rep=0;rep<n;rep++){
            ////SAVE
            H[rep]->SetTemperature(1./beta[rep]);
            PRE79MeanProperties mprop(*prop[rep],*H[rep]);
            Settings cur = settings;
            cur.hamiltonian.temperature = H[rep]->GetTemperature();//1./beta[index[i]];
            database.StoreFinalLattice(cur,lattices[rep]);
            database.StoreThermalizationHistory(cur,*prop[rep]);
            database.StoreProperties(cur,mprop,nc-1);
        }
        
        nc = settings.simulation.production_cycles;
        
        for(int rep=0;rep<n;rep++){
            prop.pop_back();
            ac.pop_back();
            simulations.pop_back();
        }
        for(int rep=0;rep<n;rep++){
            prop.push_back(new PRE79StandardProperties(&lattices[rep],nc/settings.simulation.measure_frequency));
            ac.push_back(new AutoCorrelationTimeCalculator(&lattices[rep],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length)); 
            simulations.push_back(new LatticeSimulation(H[rep],&lattices[rep],m[rep],nc));
        }
        
        Log() << "Production\n";
        
        InnerLoop(nc);
        
        for(int rep=0;rep<n;rep++){
            ////SAVE
            H[rep]->SetTemperature(1./beta[rep]);
            PRE79MeanProperties mprop(*prop[rep],*H[rep]);
            Settings cur = settings;
            cur.hamiltonian.temperature = H[rep]->GetTemperature();//1./beta[index[i]];
            database.StoreFinalLattice(cur,lattices[rep]);
            database.StorePropertiesEvolution(cur,*prop[rep]);
            database.StoreFinalProperties(cur,mprop);
        }
        
        Log() << "Finished\n";
        
        
        
        
    }
    
    void Run3() {
        int nc = settings.simulation.thermalization_cycles;
        int swapfq = settings.scanning.parallel_tempering_swapfq;
        omp_set_dynamic(0);
        omp_set_num_threads(n);
        
        for(int i=0;i<n;i++)
        {
            
            H.push_back(new PRE79StandardHamiltonian(1./beta[i],settings.hamiltonian.lambda,settings.hamiltonian.tau,settings.hamiltonian.h));
            m.push_back(new Metropolis(settings,H.back()));
            prop.push_back(new PRE79StandardProperties(&lattices[i],nc/settings.simulation.measure_frequency));
            simulations.push_back(new LatticeSimulation(H.back(),&lattices[i],m.back(),nc));
            ac.push_back(new AutoCorrelationTimeCalculator(&lattices[i],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length));
        }
        

        for(int cycle=0;cycle<nc;cycle++){
        
            //int rep = omp_get_thread_num();
            #pragma omp parallel for
            for(int rep=0;rep<n;rep++){
            //for(int cycle=0;cycle<nc;cycle++){
                ////MONTECARLO
                H[rep]->SetTemperature(1./beta[rep]);
                simulations[rep]->SetLattice(&lattices[rep]);
                prop[rep]->SetLattice(&lattices[rep]);
                
                simulations[rep]->Iterate();
                ac[rep]->Update();
                if(cycle%settings.simulation.measure_frequency==0)
                    prop[rep]->Update(cycle/settings.simulation.measure_frequency,H[rep],ac[rep]);

                        
                if(cycle==(nc-1)){
                    ////SAVE
                    H[rep]->SetTemperature(1./beta[rep]);
                    PRE79MeanProperties mprop(*prop[rep],*H[rep]);
                    Settings cur = settings;
                    cur.hamiltonian.temperature = H[rep]->GetTemperature();//1./beta[index[i]];
                    database.StoreFinalLattice(cur,lattices[rep]);
                    database.StoreThermalizationHistory(cur,*prop[rep]);
                    database.StoreProperties(cur,mprop,nc-1);
                }
            }
            
            //#pragma omp master
            if(cycle%swapfq==0){
                ////SWAP
                
                for(int rep=0;rep<n-1;rep+=2)
                    Swap(rep+1,rep,cycle/swapfq);
                for(int rep=1;rep<n-1;rep+=2)
                    Swap(rep+1,rep,cycle/swapfq);
                
                Log() << "Energies: " << E << std::endl;
                Log() << "New indices: " << index << std::endl;
                Log() << "Acceptance rate: " << acc << std::endl;
                
                ////ADJUST BETA
                
                double aptm=0.0;
                                      
                for(int i=1;i<n;i++)
                    aptm+=acc[i-1]*(beta[i]-beta[i-1]);

                double lambda = (beta[n-1]-beta[0])/aptm;

                Log() << "lambda = " << lambda << std::endl;

                vect oldbeta = beta;
                for(int i=1;i<n;i++){
                    beta[i]=beta[i-1]+lambda*acc[i-1]*(oldbeta[i]-oldbeta[i-1]);
                    /*if(beta[i]==beta[i-1]){
                        beta[i-1]*=0.95;
                        beta[i]*=1.05;
                    }*/
                }
                Log() << "New betas: " << beta << std::endl;
                
            }
            
        }
        
    }
    
    /*
    ///zla paralelizacja!?
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
                if(k%settings.simulation.measure_frequency==0){
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
                   Log() << "Simulation " << i << " has Beta " << / *beta[index[i]]* / beta[i] << " but index " << index[i] << std::endl;
                   Log() << "No of accepted swaps: " << acc << std::endl;
                   
                   std::cout << k/swapfq << " " << index << std::endl;
                }
              
                //#pragma omp barrier
                / *
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
                }* /
                 
               
                
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
                   Log() << "Simulation " << i << " has Beta " << / *beta[index[i]]* / beta[i] << " but index " << index[i] << std::endl;
                   Log() << "No of accepted swaps: " << acc << std::endl;
                   
                }
              
                //#pragma omp barrier
                / *
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
                }* /
                 
               
                
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
    }*/
private:
    void InnerLoop(int nc){
        //co ile cykli zmieniamy konfiguracje
        int swapfq = settings.scanning.parallel_tempering_swapfq;
        //replika, która teraz będzie odpowiadała za ustalanie beta
        int swaprep = 0;
        
        #pragma omp parallel for schedule(runtime) private(random01) default(shared)
        for(int rep=0;rep<n;rep++){
        
            m[rep]->AdjustRadius(&lattices[rep]);
            
            for(int cycle=0;cycle<nc;cycle++){
                ////MONTECARLO
                H[rep]->SetTemperature(1./beta[rep]);
                simulations[rep]->SetLattice(&lattices[rep]);
                prop[rep]->SetLattice(&lattices[rep]);
                
                simulations[rep]->Iterate();
                ac[rep]->Update();
                if(cycle%settings.simulation.measure_frequency==0)
                    prop[rep]->Update(cycle/settings.simulation.measure_frequency,H[rep],ac[rep]);

                        

            
            
                //#pragma omp master
                if(cycle%swapfq==0){
                    ////SWAP
                    /*
                     * Tutaj chodzi o to, aby zamieniać tylko konfiguracje sąsiadujące w beta. 
                     * Ponieważ to, czy będziemy w drugim kroku porównywali dwie sąsiadujące bądź nie,
                     * trzeba uwzględnić wiedzę, czy zamiana zaszła, czy nie. Robimy dwa kroki dlatego,
                     * żeby nie było obciążenia błądzenia przypadkowego w jedną ze stron.
                     */
                    if(rep>0 && rep<(n-1)){
                        if(Swap(rep,rep-1,cycle/swapfq/2))
                            Swap(rep+1,rep-1,cycle/swapfq/2,rep+1,rep);
                        else
                            Swap(rep+1,rep,cycle/swapfq/2);
                    }
                    else {
                        if(rep<(n-1)) Swap(rep+1,rep,cycle/swapfq);
                        if(rep>0) Swap(rep,rep-1,cycle/swapfq);
                    }
                    
                    /*for(int rep=0;rep<n-1;rep+=2)
                        Swap(rep+1,rep,cycle/swapfq);
                    for(int rep=1;rep<n-1;rep+=2)
                        Swap(rep+1,rep,cycle/swapfq);
                     */

                    Log() << "Energies: " << E << std::endl;
                    Log() << "New indices: " << index << std::endl;
                    Log() << "Acceptance rate: " << acc << std::endl;
                }
                if(cycle%swapfq==0 and cycle>3*swapfq and rep==swaprep){
                    ////ADJUST BETA
                    /*
                     * Algorytm uwzględniający średnie prawdopodobieństwa zamian na połączeniach temperatur.
                     * Takie uzgadnianie parametru beta gwarantuje podobne przekrywanie histogramów, które
                     * gwarantuje prawdziwe błądzenie przypadkowe w przestrzeni beta.
                     * Zmodyfikowana wersja algorytmu: Kerler,Rehberg, PRE 50 no. 7, pp. 4220-4225 (1994)
                     */
                    double aptm=0.0;

                    for(int i=1;i<n;i++)
                        aptm+=acc[i-1]*(beta[i]-beta[i-1]);

                    double lambda = (beta[n-1]-beta[0])/aptm;

                    Log() << "lambda = " << lambda << std::endl;

                    vect oldbeta = beta;
                    for(int i=1;i<n;i++)
                        beta[i]=beta[i-1]+lambda*acc[i-1]*(oldbeta[i]-oldbeta[i-1]);

                    Log() << "New betas: " << beta << std::endl;
                    
                    
                    //to gwarantuje, że każda replika spędzi mniej więcej tyle samo czasu na ustalaniu beta
                    if(swaprep<(n-1))
                        swaprep++;
                    else swaprep=0;
                }
            }
            
        }
    }
    
    void UpdateE(int i){
        vect e;
        int N = lattices[i].GetN();
        e.resize(N,0.0);
        for(int k=0;k<N;k++)
            e[k]=lattices[i].GetParticles()[k].GetEnergy();
        
        E[i]=e.sum()/N;
    }
    bool Swap(const int & U1, const int & U2,const int &k,const int & a1=-1, const int & a2=-1) {
        int u1 = U1;
        int u2 = U2;
        int w1 = U1;
        int w2 = U2;
        if(a1!=-1) w1=a1;
        if(a2!=-1) w2=a2;
        
        UpdateE(u1);
        UpdateE(u2);
        double dB = beta[w1]-beta[w2];
        double dE = E[u1] - E[u2];
        double frac = std::exp(dB*dE);
        if(frac>=1.0) frac=1.0;
        
        double delta = frac - acc[u2];
        acc[u2]+=delta/(k+1);

        std::cout << u1 << "->" << u2 << ", dB=" << dB << " dE=" << dE << ", exp(-dBdE)=" << frac << std::endl;

        if(random01()<frac){
            Lattice buf = lattices[u1];
            
            lattices[u1]=lattices[u2];
            #pragma omp flush(lattices)
            lattices[u2]=buf;
            #pragma omp flush(lattices)
            
            UpdateE(u1);
            UpdateE(u2);

            int ii = index[u1];
            index[u1]=index[u2];
            index[u2]=ii;
            
            return true;
             
        }else return false;
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
    std::ostream & Log() {
        ILoggable::Log() << "Thread:" << omp_get_thread_num() << "/" << omp_get_num_threads() << ": " ;
        return ILoggable::Log();
    }
};

#endif	/* PARALLELTEMPERING_H */


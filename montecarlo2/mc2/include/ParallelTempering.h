/* 
 * File:   ParallelTempering.h
 * Author: karolnew
 *
 * Created on 3 grudzień 2011, 18:13
 */

#ifndef PARALLELTEMPERING_H
#define	PARALLELTEMPERING_H

#include "ILoggable.h"
#include "Statistical.h"
#include "omp.h"


class ParallelTempering: public ILoggable {
    std::vector<shared_ptr<LatticeSimulation> >      simulations;
    std::vector<shared_ptr<Lattice> >                lattices;
    std::vector<shared_ptr<PRE79StandardProperties> > prop;
    std::vector<shared_ptr<PRE79StandardHamiltonian> > H;
    std::vector<shared_ptr<Metropolis> >     m;
    std::vector<shared_ptr<AutoCorrelationTimeCalculator> > ac;
    
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
                lattices.push_back(shared_ptr<Lattice>(new Lattice(settings,
                        Lattice::Isotropic)));
        else if(settings.initial.isotropic && settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(shared_ptr<Lattice>(new Lattice(settings,
                        Lattice::IsotropicRighthanded)));
        else if(settings.initial.biaxial && !settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(shared_ptr<Lattice>(new Lattice(settings,
                        Lattice::Biaxial)));
        else if(settings.initial.biaxial && settings.initial.righthanded)
        for(int i=0;i<n;i++)
                lattices.push_back(shared_ptr<Lattice>(new Lattice(settings,
                        Lattice::BiaxialRighthanded)));
        
            
            
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
            
            H.push_back(shared_ptr<PRE79StandardHamiltonian>(new PRE79StandardHamiltonian(1./beta[i],settings.hamiltonian.lambda,settings.hamiltonian.tau,settings.hamiltonian.h)));
            m.push_back(shared_ptr<Metropolis>(new Metropolis(settings,dynamic_pointer_cast<Hamiltonian>(H[i]))));
            prop.push_back(shared_ptr<PRE79StandardProperties>(new PRE79StandardProperties(lattices[i],nc/settings.simulation.measure_frequency)));
            simulations.push_back(shared_ptr<LatticeSimulation>(new LatticeSimulation(dynamic_pointer_cast<Hamiltonian>(H[i]),lattices[i],m[i],nc)));
            ac.push_back(shared_ptr<AutoCorrelationTimeCalculator>(new AutoCorrelationTimeCalculator(lattices[i],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length)));
        }
        Log() << "Thermalization\n";
        InnerLoop(nc);
        
        Log() << "Saving thermalization properties\n";
        
        
        /*
         * "prop" zawierają teraz mieszankę właściwości systemu dla różnych temperatur
         * i tak naprawdę dla różnych replik. Trzeba rozdysponować te wartości
         * do poprawnych tablic.
         */
        
        for(int rep=0;rep<n;rep++){
            shared_ptr<PRE79StandardProperties> proper_prop(new PRE79StandardProperties(lattices[rep],0));
            foreach(shared_ptr<PRE79StandardProperties> p,prop){
                proper_prop->Append(PRE79StandardProperties::FromReplicaMask(*p,(double)rep));
            }
            
            H[rep]->SetTemperature(1./beta[rep]);
            PRE79MeanProperties mprop(proper_prop,H[rep]);
            Settings cur = settings;
            cur.hamiltonian.temperature = H[rep]->GetTemperature();//1./beta[index[i]];
            database.StoreFinalLattice(cur,*lattices[rep]);
            database.StoreThermalizationHistory(cur,*proper_prop);
            database.StoreProperties(cur,mprop,nc-1);
        }

        
        nc = settings.simulation.production_cycles;
        
        prop.clear();
        ac.clear();
        simulations.clear();

        for(int rep=0;rep<n;rep++){
            prop.push_back(shared_ptr<PRE79StandardProperties>(new PRE79StandardProperties(lattices[rep],nc/settings.simulation.measure_frequency)));
            ac.push_back(shared_ptr<AutoCorrelationTimeCalculator>(new AutoCorrelationTimeCalculator(lattices[rep],settings.simulation.autocorrelation_frequency,settings.simulation.autocorrelation_length))); 
            simulations.push_back(shared_ptr<LatticeSimulation>(new LatticeSimulation(dynamic_pointer_cast<Hamiltonian>(H[rep]),lattices[rep],m[rep],nc)));
        }
        
        Log() << "Production\n";
        
        InnerLoop(nc);
        
        for(int rep=0;rep<n;rep++){
            shared_ptr<PRE79StandardProperties> proper_prop(new PRE79StandardProperties(lattices[rep],0));
            foreach(shared_ptr<PRE79StandardProperties> p,prop){
                proper_prop->Append(PRE79StandardProperties::FromReplicaMask(*p,(double)rep));
            }
            
            proper_prop->CalculateSpecificHeat();
            
            H[rep]->SetTemperature(1./beta[rep]);
            PRE79MeanProperties mprop(proper_prop,H[rep]);
            Settings cur = settings;
            cur.hamiltonian.temperature = H[rep]->GetTemperature();//1./beta[index[i]];
            database.StoreFinalLattice(cur,*lattices[rep]);
            database.StorePropertiesEvolution(cur,*proper_prop);
            database.StoreFinalProperties(cur,mprop);
        }

        
        Log() << "Finished\n";
        
        
        
        
    }
   
private:
    void InnerLoop(int nc){
        //co ile cykli zmieniamy konfiguracje
        int swapfq = settings.scanning.parallel_tempering_swapfq;
        //replika, która teraz będzie odpowiadała za ustalanie beta
        //int swaprep = 0;
        
        //#pragma omp parallel for schedule(runtime) private(random01)
        //for(int rep=0;rep<n;rep++){
        
        #pragma omp parallel private(random01)
        {
            int rep = omp_get_thread_num();
        
            m[rep]->AdjustRadius(lattices[rep]);
            
            for(int cycle=0;cycle<nc;cycle++){
                ////MONTECARLO
                //H[rep]->SetTemperature(1./beta[index[rep]]);
                //H[rep]->SetTemperature(1./prop[rep]->GetCurrentBetaValue());
                //simulations[rep]->SetLattice(lattices[rep]);
                //prop[rep]->SetLattice(lattices[rep]);
                
                simulations[rep]->Iterate();
                ac[rep]->Update();
                if(cycle%settings.simulation.measure_frequency==0)
                    prop[rep]->Update(cycle/settings.simulation.measure_frequency,H[rep],ac[rep]);

                
                #pragma omp barrier
                
                if(cycle%swapfq==0 && cycle>0)
                #pragma omp critical
                {
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
                
//                #pragma omp master
//                if(cycle%swapfq==0 and cycle>3*swapfq){
//                    ////ADJUST BETA
//                    /*
//                     * Algorytm uwzględniający średnie prawdopodobieństwa zamian na połączeniach temperatur.
//                     * Takie uzgadnianie parametru beta gwarantuje podobne przekrywanie histogramów, które
//                     * gwarantuje prawdziwe błądzenie przypadkowe w przestrzeni beta.
//                     * Zmodyfikowana wersja algorytmu: Kerler,Rehberg, PRE 50 no. 7, pp. 4220-4225 (1994)
//                     */
//                    double aptm=0.0;
//
//                    for(int i=1;i<n;i++)
//                        aptm+=acc[i-1]*(beta[i]-beta[i-1]);
//
//                    double lambda = (beta[n-1]-beta[0])/aptm;
//
//                    Log() << "lambda = " << lambda << std::endl;
//
//                    vect oldbeta = beta;
//                    for(int i=1;i<n;i++)
//                        beta[i]=beta[i-1]+lambda*acc[i-1]*(oldbeta[i]-oldbeta[i-1]);
//
//                    Log() << "New betas: " << beta << std::endl;
//                    
//                    
//                    /*
//                    //to gwarantuje, że każda replika spędzi mniej więcej tyle samo czasu na ustalaniu beta
//                    if(swaprep<(n-1))
//                        swaprep++;
//                    else swaprep=0;
//                     */
//                }
            }
            
        }
    }
    
    void UpdateE(int i){
        vect e;
        int N = lattices[i]->GetN();
        e.resize(N,0.0);
        for(int k=0;k<N;k++)
            e[k]=lattices[i]->GetParticles()[k].GetEnergy();
        
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

        Log() << u1 << "->" << u2 << ", dB=" << dB << " dE=" << dE << ", exp(-dBdE)=" << frac << std::endl;

        if(random01()<frac){
            Log() << "Swapping!\n";

            
            prop[u1]->SetCurrentBetaValue(beta[index[u2]]);
            prop[u1]->SetCurrentBetaIndex(index[u2]);
            prop[u2]->SetCurrentBetaValue(beta[index[u1]]);
            prop[u2]->SetCurrentBetaIndex(index[u1]);
            
            H[u1]->SetTemperature(1./beta[index[u2]]);
            H[u2]->SetTemperature(1./beta[index[u1]]);
            
            
            UpdateE(u1);
            UpdateE(u2);

            int ii = index[u1];
            index[u1]=index[u2];
            index[u2]=ii;
            
            return true;
             
        }else return false;
    }
    
    std::ostream & Log() {
        ILoggable::Log() << "Thread:" << omp_get_thread_num() << "/" << omp_get_num_threads() << ": " ;
        return ILoggable::Log();
    }
};

#endif	/* PARALLELTEMPERING_H */


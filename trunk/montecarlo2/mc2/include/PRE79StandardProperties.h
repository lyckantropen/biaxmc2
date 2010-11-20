/* 
 * File:   StandardProperties.h
 * Author: karol
 *
 * Created on 18 listopad 2009, 11:27
 */

#ifndef _STANDARDPROPERTIES_H
#define	_STANDARDPROPERTIES_H

//#include "SpatialCorrelation.h"
#include "SpatialCorrelationEvolution.h"
#include "PRE79StandardHamiltonian.h"
#include "Lattice.h"
#include "Statistical.h"
#include "Contractions.h"
#include "PRE79SpatialCorrelations.h"

#include "serializer.h"
#include "eig3.h"
#include "evsort.h"


//extern "C" void Jacobi_Cyclic_Method(double *eigenvalues, double *eigenvectors, double *A, int n);


///obliczanie właściwości układu i przechowywanie historii stanów
class PRE79StandardProperties {
    template<class serializer_t>
    friend void operator|(serializer_t & s, PRE79StandardProperties & prop);

    Lattice *   lat;        ///<wskaźnik do bieżącej siatki
    bool    readonly;       ///<flaga sygnalizująca tryb tylko do odczytu
    int index;              ///<oznaczenie właściwości - numer cyklu produkcyjnego MC
    int acc_idx;            ///<numer chwili czasu
    int ncycles;            ///<liczba cykli do średniej
    
    vect    energy;         ///<energia w funkcji czasu
    vect    parity;         ///<średnia parzystość w funkcji czasu
//    vect    D200;
//    vect    D220;
//    vect    D202;
//    vect    D222;
    //Value   specific_heat;  ///<ciepło właściwe w funkcji czasu (w praktyce obliczanie całej historii jest bardzo kosztowne, obliczamy dopiero na końcu)
    Value   fluctuation;    ///<fluktuacje energii


    //ewolucja średnich tensorów xx,yy,zz
    std::vector<vect>  MeanQxTensor;
    std::vector<vect>  MeanQyTensor;
    std::vector<vect>  MeanQzTensor;

    //ewolucja funkcji korelacji dla różnych czasów
    Delta200CorrelationZ d200corz;
    Delta222CorrelationZ d222corz;
    Delta220CorrelationZ d220corz;

    Delta200CorrelationX d200corx;
    Delta222CorrelationX d222corx;
    Delta220CorrelationX d220corx;

    Delta200CorrelationY d200cory;
    Delta222CorrelationY d222cory;
    Delta220CorrelationY d220cory;

    Delta322Correlation d322cor;
    ParityCorrelation paritycor;

    ///obliczanie średniej energii (po siatce)
    void CalculateMeanEnergy(){
        double E=0.0;
        for(int site=0;site<lat->GetN();site++){
            E+=lat->GetParticles()[site].GetEnergy();
        }
        E/=double(lat->GetN());
        energy[acc_idx]=E;
    }
    ///obliczanie średniej parzystości
    void CalculateMeanParity(){
        double p=0.0;
        for(int site=0;site<lat->GetN();site++)
            p+=lat->GetParticles()[site].GetParity();
        p/=double(lat->GetN());
        parity[acc_idx]=p;
    }
public:
    ///obliczanie ciepła właściwego metodą bootstrapu
    void CalculateSpecificHeat(){
        fluctuation=CalculateFluctuation(energy,acc_idx)*Value(lat->GetN());
    }
private:
    ///obliczanie średnich tensorów xx,yy i zz
    void CalculateMeanTensors(){
        vect mqx(0.0,6);
        vect mqy(0.0,6);
        vect mqz(0.0,6);
        for(int i=0;i<lat->GetN();i++){
            mqx+=lat->GetParticles()[i].GetQX();
            mqy+=lat->GetParticles()[i].GetQY();
            mqz+=lat->GetParticles()[i].GetQZ();
        }
        MeanQxTensor[acc_idx]=mqx/double(lat->GetN());
        MeanQyTensor[acc_idx]=mqy/double(lat->GetN());
        MeanQzTensor[acc_idx]=mqz/double(lat->GetN());
    }

public:
    ///konstruktor serializacyjny
    PRE79StandardProperties():readonly(true),index(0),acc_idx(0) {}
    ///konstruktor tradycyjny
    PRE79StandardProperties(Lattice * l, int _ncycles):
    lat(l),readonly(false),ncycles(_ncycles),
    d200corz(l,_ncycles),
    d222corz(l,_ncycles),
    d220corz(l,_ncycles),
    d200corx(l,_ncycles),
    d222corx(l,_ncycles),
    d220corx(l,_ncycles),
    d200cory(l,_ncycles),
    d222cory(l,_ncycles),
    d220cory(l,_ncycles),
    d322cor(l,_ncycles),
    MeanQxTensor(_ncycles),
    MeanQyTensor(_ncycles),
    MeanQzTensor(_ncycles),
//    D200(0.0,_ncycles),
//    D220(0.0,_ncycles),
//    D202(0.0,_ncycles),
//    D222(0.0,_ncycles),
    paritycor(l,_ncycles)
    {
        acc_idx=-1;
        index=0;
        energy.resize(ncycles,0.0);
        parity.resize(ncycles,0.0);
        for(int i=0;i<ncycles;i++){
            MeanQxTensor[i].resize(6,0.0);
            MeanQyTensor[i].resize(6,0.0);
            MeanQzTensor[i].resize(6,0.0);
        }
        //specific_heat.resize(ncycles);
    }


    PRE79StandardProperties(const PRE79StandardProperties & p){
        //tak to musi być zrobione, por. Standard C++ Library Reference pp. 327
//        D200.resize(p.D200.size());
//        D220.resize(p.D220.size());
//        D202.resize(p.D202.size());
//        D222.resize(p.D222.size());
        energy.resize(p.energy.size());
        energy = p.energy;
        parity.resize(p.parity.size());
        parity = p.parity;
        lat=p.lat;
        readonly=p.readonly;
        index=p.index;
        ncycles=p.ncycles;
        //specific_heat=p.specific_heat;
        fluctuation=p.fluctuation;
        d200corz=p.d200corz;
        d220corz=p.d220corz;
        d222corz=p.d222corz;
        d200corx=p.d200corx;
        d220corx=p.d220corx;
        d222corx=p.d222corx;
        d200cory=p.d200cory;
        d220cory=p.d220cory;
        d222cory=p.d222cory;
        MeanQxTensor=p.MeanQxTensor;
        MeanQyTensor=p.MeanQyTensor;
        MeanQzTensor=p.MeanQzTensor;
        d322cor=p.d322cor;
        paritycor=p.paritycor;
        acc_idx=p.acc_idx;
    }
    const PRE79StandardProperties operator=(const PRE79StandardProperties & p){
        /*
        D200.resize(p.D200.size());
        D220.resize(p.D220.size());
        D202.resize(p.D202.size());
        D222.resize(p.D222.size());
        */
        energy.resize(p.energy.size());
        energy = p.energy;
        parity.resize(p.parity.size());
        parity = p.parity;
        lat=p.lat;
        readonly=p.readonly;
        index=p.index;
        ncycles=p.ncycles;
        //specific_heat=p.specific_heat;
        fluctuation=p.fluctuation;
        d200corz=p.d200corz;
        d220corz=p.d220corz;
        d222corz=p.d222corz;
        d200corx=p.d200corx;
        d220corx=p.d220corx;
        d222corx=p.d222corx;
        d200cory=p.d200cory;
        d220cory=p.d220cory;
        d222cory=p.d222cory;
        MeanQxTensor=p.MeanQxTensor;
        MeanQyTensor=p.MeanQyTensor;
        MeanQzTensor=p.MeanQzTensor;
        d322cor=p.d322cor;
        paritycor=p.paritycor;
        acc_idx=p.acc_idx;
        return *this;
    }

    ///obliczenie kontrakcji
    void Update(int idx,const PRE79StandardHamiltonian * H){
        if(readonly) return;
        if((acc_idx+1)>=ncycles) return;
        index=idx;
        acc_idx++;

        d200corz.Update();
        d220corz.Update();
        d222corz.Update();
        d200corx.Update();
        d220corx.Update();
        d222corx.Update();
        d200cory.Update();
        d220cory.Update();
        d222cory.Update();

        d322cor.Update();
        paritycor.Update();

        CalculateMeanEnergy();
        CalculateMeanParity();
        CalculateMeanTensors();

        if((acc_idx+1)>=ncycles)
            CalculateSpecificHeat();
    }

    ///dołożenie danych z innych właściwości, koniecznie po zakończeniu obliczeń
    void Append(const PRE79StandardProperties & p){
        ncycles+=p.ncycles;
        acc_idx+=p.ncycles;
        index+=p.ncycles;
        d200corz.Append(p.d200corz);
        d220corz.Append(p.d220corz);
        d222corz.Append(p.d222corz);
        d200corx.Append(p.d200corx);
        d220corx.Append(p.d220corx);
        d222corx.Append(p.d222corx);
        d200cory.Append(p.d200cory);
        d220cory.Append(p.d220cory);
        d222cory.Append(p.d222cory);
        d322cor.Append(p.d322cor);
        paritycor.Append(p.paritycor);

        MeanQxTensor.insert(MeanQxTensor.end(),p.MeanQxTensor.begin(),p.MeanQxTensor.end());
        MeanQyTensor.insert(MeanQyTensor.end(),p.MeanQyTensor.begin(),p.MeanQyTensor.end());
        MeanQzTensor.insert(MeanQzTensor.end(),p.MeanQzTensor.begin(),p.MeanQzTensor.end());

        int oldsize=energy.size();
        vect tmpenergy(0.0,ncycles);
        vect tmpparity(0.0,ncycles);
        /*
        vect tmpD200(0.0,ncycles);
        vect tmpD220(0.0,ncycles);
        vect tmpD202(0.0,ncycles);
        vect tmpD222(0.0,ncycles);
         */
        for(int i=0;i<oldsize;i++){
            tmpenergy[i]=energy[i];
            tmpparity[i]=parity[i];
            /*
            tmpD200[i]=D200[i];
            tmpD220[i]=D220[i];
            tmpD202[i]=D202[i];
            tmpD222[i]=D222[i];
             */
        }

        for(int i=oldsize;i<ncycles;i++){
            tmpenergy[i]=p.energy[i-oldsize];
            tmpparity[i]=p.parity[i-oldsize];
            /*
            tmpD200[i]=p.D200[i-oldsize];
            tmpD220[i]=p.D220[i-oldsize];
            tmpD202[i]=p.D202[i-oldsize];
            tmpD222[i]=p.D222[i-oldsize];
             */
        }
        energy.resize(ncycles,0.0);
        parity.resize(ncycles,0.0);
        //D200.resize(ncycles,0.0);
        //D220.resize(ncycles,0.0);
        //D202.resize(ncycles,0.0);
        //D222.resize(ncycles,0.0);
        energy=tmpenergy;
        parity=tmpparity;
        //D200=tmpD200;
        //D220=tmpD220;
        //D202=tmpD202;
        //D222=tmpD222;
    }

    Value TemporalMeanEnergyPerMolecule() const {
        //return BootstrapMean(energy,0,acc_idx+1);
        return Mean(energy,0,acc_idx+1);
    }
    Value TemporalMeanParity() const {
        return Mean(parity,0,acc_idx+1);
    }
    Value ParitySusceptibility() const {
        return CalculateFluctuation(parity,acc_idx)*Value(lat->GetN());
    }
    /*
    Value MeanDelta200() const {
        return BootstrapMean(D200,0,acc_idx+1);
    }
    Value MeanDelta220() const {
        return BootstrapMean(D220,0,acc_idx+1);
    }
    Value MeanDelta202() const {
        return BootstrapMean(D202,0,acc_idx+1);
    }
    Value MeanDelta222() const {
        return BootstrapMean(D222,0,acc_idx+1);
    }
     */

    //korelacje dla różnych osi
    Value Delta200ZByCorrelation() const {
        return sqrt(d200corz.Limit());
    }
    Value Delta200ZByCorrelationSusceptibility() const {
        return CalculateFluctuation(d200corz.LimitHistory(),acc_idx)*Value(lat->GetN());
    }
    Value Delta222ZByCorrelation() const {
        return sqrt(d222corz.Limit());
    }
    Value Delta222ZByCorrelationSusceptibility() const {
        return CalculateFluctuation(d222corz.LimitHistory(),acc_idx)*Value(lat->GetN());
    }
    Value Delta200XByCorrelation() const {
        return sqrt(d200corx.Limit());
    }
    Value Delta200XByCorrelationSusceptibility() const {
        return CalculateFluctuation(d200corx.LimitHistory(),acc_idx)*Value(lat->GetN());
    }
    Value Delta222XByCorrelation() const {
        return sqrt(d222corx.Limit());
    }
    Value Delta222XByCorrelationSusceptibility() const {
        return CalculateFluctuation(d222corx.LimitHistory(),acc_idx)*Value(lat->GetN());
    }
    Value Delta200YByCorrelation() const {
        return sqrt(d200cory.Limit());
    }
    Value Delta200YByCorrelationSusceptibility() const {
        return CalculateFluctuation(d200cory.LimitHistory(),acc_idx)*Value(lat->GetN());
    }
    Value Delta222YByCorrelation() const {
        return sqrt(d222cory.Limit());
    }
    Value Delta222YByCorrelationSusceptibility() const {
        return CalculateFluctuation(d222cory.LimitHistory(),acc_idx)*Value(lat->GetN());
    }

    Value Delta322ByCorrelation() const {
        return sqrt(d322cor.Limit());
    }
    Value Delta322ByCorrelationSusceptibility() const {
        return CalculateFluctuation(d322cor.LimitHistory(),acc_idx)*Value(lat->GetN());
    }
    Value ParityByCorrelation() const {
        return sqrt(paritycor.Limit());
    }
    Value ParityByCorrelationSusceptibility() const {
        return CalculateFluctuation(paritycor.LimitHistory(),acc_idx)*Value(lat->GetN());
    }

    //średnie funkcje korelacji dla poszczególnych osi
    vect Delta200ZMeanCorrelation() const {
//        std::cout << "D200Z\n";
        return d200corz.Mean(acc_idx+1);
    }
    vect Delta222ZMeanCorrelation() const {
//        std::cout << "D222Z\n";
        return d222corz.Mean(acc_idx+1);
    }
    vect Delta220ZMeanCorrelation() const {
//        std::cout << "D220Z\n";
        return d220corz.Mean(acc_idx+1);
    }
    vect Delta200YMeanCorrelation() const {
//        std::cout << "D200Y\n";
        return d200cory.Mean(acc_idx+1);
    }
    vect Delta222YMeanCorrelation() const {
//        std::cout << "D222Y\n";
        return d222cory.Mean(acc_idx+1);
    }
    vect Delta220YMeanCorrelation() const {
//        std::cout << "D220Y\n";
        return d220cory.Mean(acc_idx+1);
    }
    vect Delta200XMeanCorrelation() const {
//        std::cout << "D200X\n";
        return d200corx.Mean(acc_idx+1);
    }
    vect Delta222XMeanCorrelation() const {
//        std::cout << "D222X\n";
        return d222corx.Mean(acc_idx+1);
    }
    vect Delta220XMeanCorrelation() const {
//        std::cout << "D220X\n";
        return d220corx.Mean(acc_idx+1);
    }



    vect Delta322MeanCorrelation() const {
//        std::cout << "D322\n";
        return d322cor.Mean(acc_idx+1);
    }
    vect ParityMeanCorrelation() const {
//        std::cout << "Parity\n";
        return paritycor.Mean(acc_idx+1);
    }


    const vect & GetDelta200ZCorrelation() const {
        return d200corz.CurrentCorrelation();
    }
    const vect & GetDelta222ZCorrelation() const {
        return d222corz.CurrentCorrelation();
    }
    const vect & GetDelta220ZCorrelation() const {
        return d220corz.CurrentCorrelation();
    }
    const vect & GetDelta200XCorrelation() const {
        return d200corx.CurrentCorrelation();
    }
    const vect & GetDelta222XCorrelation() const {
        return d222corx.CurrentCorrelation();
    }
    const vect & GetDelta220XCorrelation() const {
        return d220corx.CurrentCorrelation();
    }
    const vect & GetDelta200YCorrelation() const {
        return d200cory.CurrentCorrelation();
    }
    const vect & GetDelta222YCorrelation() const {
        return d222cory.CurrentCorrelation();
    }
    const vect & GetDelta220YCorrelation() const {
        return d220cory.CurrentCorrelation();
    }


    const vect & GetDelta322Correlation() const {
        return d322cor.CurrentCorrelation();
    }
    const vect & GetParityCorrelation() const {
        return paritycor.CurrentCorrelation();
    }

    const vect GetMeanQxTensor() const {
//        std::cout << "Qx\n";
        return MeanVector(MeanQxTensor,0,acc_idx+1);
    }
    const vect GetMeanQyTensor() const {
//        std::cout << "Qy\n";
        return MeanVector(MeanQyTensor,0,acc_idx+1);
    }
    const vect GetMeanQzTensor() const {
//        std::cout << "Qz\n";
        return MeanVector(MeanQzTensor,0,acc_idx+1);
    }
    
    const vect & EnergyEvolution() const {
        return energy;
    }
    const vect & ParityEvolution() const {
        return parity;
    }

    const Delta200CorrelationZ & Delta200ZCorrelationEvolution() const {
        return d200corz;
    }
    const Delta222CorrelationZ & Delta222ZCorrelationEvolution() const {
        return d222corz;
    }
    const Delta220CorrelationZ & Delta220ZCorrelationEvolution() const {
        return d220corz;
    }
    const Delta200CorrelationX & Delta200XCorrelationEvolution() const {
        return d200corx;
    }
    const Delta222CorrelationX & Delta222XCorrelationEvolution() const {
        return d222corx;
    }
    const Delta220CorrelationX & Delta220XCorrelationEvolution() const {
        return d220corx;
    }
    const Delta200CorrelationY & Delta200YCorrelationEvolution() const {
        return d200cory;
    }
    const Delta222CorrelationY & Delta222YCorrelationEvolution() const {
        return d222cory;
    }
    const Delta220CorrelationY & Delta220YCorrelationEvolution() const {
        return d220cory;
    }

    const Delta322Correlation & Delta322CorrelationEvolution() const {
        return d322cor;
    }
    
    const ParityCorrelation & ParityCorrelationEvolution() const {
        return paritycor;
    }

    //Accessors
    ///ciepło właściwe
    const Value & SpecificHeat() const {
        //CalculateSpecificHeat(H->GetTemperature());
        //return specific_heat;
        return fluctuation;
    }
 
    const Value & Fluctuation() const {
        return fluctuation;
    }
    const int & GetIndex() const {
        return index;
    }
    const int & GetAccIdx() const {
        return acc_idx;
    }
    const int & GetNCycles() const {
        return ncycles;
    }
    int GetMaxCorrLen() const {
        return lat->GetL()/2+1;
    }

};

template<class serializer_t>
void operator|(serializer_t & s, PRE79StandardProperties & prop){
    s|prop.index;
    s|prop.acc_idx;
    s|prop.ncycles;
    //s|prop.specific_heat;
    s|prop.fluctuation;
    s|prop.d200corz;
    s|prop.d220corz;
    s|prop.d222corz;
    s|prop.d200corx;
    s|prop.d220corx;
    s|prop.d222corx;
    s|prop.d200cory;
    s|prop.d220cory;
    s|prop.d222cory;
//    s|prop.D200;
//    s|prop.D220;
//    s|prop.D202;
//    s|prop.D222;
    s|prop.d322cor;
    s|prop.paritycor;
    s|prop.MeanQxTensor;
    s|prop.MeanQyTensor;
    s|prop.MeanQzTensor;
    s|prop.energy;
    s|prop.parity;
}

/**
 * używane do przechowania końcowych, uśrednionych parametrów układu
 */
class PRE79MeanProperties {
    template<class serializer_t>
    friend void operator|(serializer_t & s, PRE79MeanProperties & prop);
    Value energy;
    Value parity;
    Value parity_sus;
    //Value specific_heat;
    Value fluctuation;
    Value d200z_from_correlation;
    Value d222z_from_correlation;
    Value d200x_from_correlation;
    Value d222x_from_correlation;
    Value d200y_from_correlation;
    Value d222y_from_correlation;
    Value d322_from_correlation;
    Value parity_from_correlation;
    Value d200z_from_correlation_sus;
    Value d222z_from_correlation_sus;
    Value d200x_from_correlation_sus;
    Value d222x_from_correlation_sus;
    Value d200y_from_correlation_sus;
    Value d222y_from_correlation_sus;
    Value d322_from_correlation_sus;
    Value parity_from_correlation_sus;
    double mean_d200;
    double mean_d220;
    double mean_d202;
    double mean_d222;

    vect mean_d200corz;
    vect mean_d220corz;
    vect mean_d222corz;
    vect mean_d200corx;
    vect mean_d220corx;
    vect mean_d222corx;
    vect mean_d200cory;
    vect mean_d220cory;
    vect mean_d222cory;

    vect mean_qx;
    vect mean_qy;
    vect mean_qz;

    vect mean_d322cor;
    vect mean_paritycor;

    double temperature;
    double tau;
    double lambda;
    double h;



public:
    ///konstruktor serializacyjny
    PRE79MeanProperties(){}
    ///konstruktor tradycyjny
    PRE79MeanProperties(PRE79StandardProperties & prop, PRE79StandardHamiltonian & H){
        mean_d200corz.resize(prop.GetMaxCorrLen(),0.0);
        mean_d222corz.resize(prop.GetMaxCorrLen(),0.0);
        mean_d220corz.resize(prop.GetMaxCorrLen(),0.0);
        mean_d200corx.resize(prop.GetMaxCorrLen(),0.0);
        mean_d222corx.resize(prop.GetMaxCorrLen(),0.0);
        mean_d220corx.resize(prop.GetMaxCorrLen(),0.0);
        mean_d200cory.resize(prop.GetMaxCorrLen(),0.0);
        mean_d222cory.resize(prop.GetMaxCorrLen(),0.0);
        mean_d220cory.resize(prop.GetMaxCorrLen(),0.0);
        mean_qx.resize(6,0.0);
        mean_qz.resize(6,0.0);
        mean_qy.resize(6,0.0);

        mean_d322cor.resize(prop.GetMaxCorrLen(),0.0);
        mean_paritycor.resize(prop.GetMaxCorrLen(),0.0);


        energy = prop.TemporalMeanEnergyPerMolecule();
        parity = prop.TemporalMeanParity();
        parity_sus = prop.ParitySusceptibility();
        //specific_heat = prop.SpecificHeat();
        fluctuation = prop.Fluctuation();

        d200z_from_correlation = prop.Delta200ZByCorrelation();
        d222z_from_correlation = prop.Delta222ZByCorrelation();
        d200x_from_correlation = prop.Delta200XByCorrelation();
        d222x_from_correlation = prop.Delta222XByCorrelation();
        d200y_from_correlation = prop.Delta200YByCorrelation();
        d222y_from_correlation = prop.Delta222YByCorrelation();
        
        d200z_from_correlation_sus = prop.Delta200ZByCorrelationSusceptibility();
        d222z_from_correlation_sus = prop.Delta222ZByCorrelationSusceptibility();
        d200x_from_correlation_sus = prop.Delta200XByCorrelationSusceptibility();
        d222x_from_correlation_sus = prop.Delta222XByCorrelationSusceptibility();
        d200y_from_correlation_sus = prop.Delta200YByCorrelationSusceptibility();
        d222y_from_correlation_sus = prop.Delta222YByCorrelationSusceptibility();

        d322_from_correlation = prop.Delta322ByCorrelation();
        parity_from_correlation = prop.ParityByCorrelation();

        d322_from_correlation_sus = prop.Delta322ByCorrelationSusceptibility();
        parity_from_correlation_sus = prop.ParityByCorrelationSusceptibility();

        mean_d200corz = prop.Delta200ZMeanCorrelation();
        mean_d220corz = prop.Delta220ZMeanCorrelation();
        mean_d222corz = prop.Delta222ZMeanCorrelation();
        mean_d200corx = prop.Delta200XMeanCorrelation();
        mean_d220corx = prop.Delta220XMeanCorrelation();
        mean_d222corx = prop.Delta222XMeanCorrelation();
        mean_d200cory = prop.Delta200YMeanCorrelation();
        mean_d220cory = prop.Delta220YMeanCorrelation();
        mean_d222cory = prop.Delta222YMeanCorrelation();
        //mean_d200 = prop.MeanDelta200();
        //mean_d220 = prop.MeanDelta220();
        //mean_d202 = prop.MeanDelta202();
        //mean_d222 = prop.MeanDelta222();


        mean_d322cor = prop.Delta322MeanCorrelation();
        mean_paritycor = prop.ParityMeanCorrelation();

        mean_qx = prop.GetMeanQxTensor();
        mean_qy = prop.GetMeanQyTensor();
        mean_qz = prop.GetMeanQzTensor();

        temperature = H.GetTemperature();
        tau = H.GetTau();
        lambda = H.GetLambda();
        h = H.GetH();
    }
    //do wykonania po wczytaniu
    void CalculateMeanTensors() {

        vect mqx(0.0,6);
        vect mqy(0.0,6);
        vect mqz(0.0,6);
        mqz = MeanQxTensor();
        mqx = MeanQyTensor();
        mqy = MeanQzTensor();

	//tensory w bazie cząsteczki
        vect t20(0.0,6);
        vect t22(0.0,6);
        t20 = std::sqrt(1.5)*(mqz - Identity(3)/3.);
        t22 = (mqx-mqy)/std::sqrt(2.);

        double t20evectm[3][3];
        double t20evalm[3];
        double t22evectm[3][3];
        double t22evalm[3];

        double t20mat[3][3];
        double t22mat[3][3];

        t20mat[0][0]=t20[0];    t20mat[0][1]=t20mat[1][0]=t20[1];   t20mat[0][2]=t20mat[2][0]=t20[2];
                                t20mat[1][1]=t20[3];                t20mat[1][2]=t20mat[2][1]=t20[4];
                                                                    t20mat[2][2]=t20[5];

        t22mat[0][0]=t22[0];    t22mat[0][1]=t22mat[1][0]=t22[1];   t22mat[0][2]=t22mat[2][0]=t22[2];
                                t22mat[1][1]=t22[3];                t22mat[1][2]=t22mat[2][1]=t22[4];
                                                                    t22mat[2][2]=t22[5];

        //diagonalizacja
        //dsyevj3(t20mat,t20evectm,t20evalm);
        //dsyevj3(t22mat,t22evectm,t22evalm);
	//Jacobi_Cyclic_Method((double*)t20evalm,(double*)t20evectm,(double*)t20mat,3);
	//Jacobi_Cyclic_Method((double*)t22evalm,(double*)t22evectm,(double*)t22mat,3);
      
	eigen_decomposition(t20mat,t20evectm,t20evalm);
        eigen_decomposition(t22mat,t22evectm,t22evalm);
	
	
	std::vector<evs> t20es;
	std::vector<evs> t22es;
	t20es.push_back(evs(t20evalm[0],t20evectm[0]));
	t20es.push_back(evs(t20evalm[1],t20evectm[1]));
	t20es.push_back(evs(t20evalm[2],t20evectm[2]));
	t22es.push_back(evs(t22evalm[0],t22evectm[0]));
	t22es.push_back(evs(t22evalm[1],t22evectm[1]));
	t22es.push_back(evs(t22evalm[2],t22evectm[2]));
	std::sort(t20es.begin(),t20es.end());
	std::sort(t22es.begin(),t22es.end());
	
	//std::cout << "ev: " << t20es[0].e << " " << t20es[1].e << " " <<  t20es[2].e << std::endl;
	
        vect t20a(0.0,3);
        vect t20b(0.0,3);
        vect t20c(0.0,3);
        vect t22a(0.0,3);
        vect t22b(0.0,3);
        vect t22c(0.0,3);

	for(int i=0;i<3;i++){
		t20c[i] = t20es[0].v[i];
		t20a[i] = t20es[1].v[i];
		t20b[i] = t20es[2].v[i];
		t22c[i] = t22es[0].v[i];
		t22a[i] = t22es[1].v[i];
		t22b[i] = t22es[2].v[i];
	}


        //iloczyny tensorowe
        vect t20aa(0.0,6);
        vect t20bb(0.0,6);
        vect t20cc(0.0,6);
        vect t22aa(0.0,6);
        vect t22bb(0.0,6);
        vect t22cc(0.0,6);

        t20aa[0]=t20a[0]*t20a[0];   t20aa[1]=t20a[1]*t20a[0];   t20aa[2]=t20a[2]*t20a[0];
                                    t20aa[3]=t20a[1]*t20a[1];   t20aa[4]=t20a[2]*t20a[1];
                                                                t20aa[5]=t20a[2]*t20a[2];

        t20bb[0]=t20b[0]*t20b[0];   t20bb[1]=t20b[1]*t20b[0];   t20bb[2]=t20b[2]*t20b[0];
                                    t20bb[3]=t20b[1]*t20b[1];   t20bb[4]=t20b[2]*t20b[1];
                                                                t20bb[5]=t20b[2]*t20b[2];

        t20cc[0]=t20c[0]*t20c[0];   t20cc[1]=t20c[1]*t20c[0];   t20cc[2]=t20c[2]*t20c[0];
                                    t20cc[3]=t20c[1]*t20c[1];   t20cc[4]=t20c[2]*t20c[1];
                                                                t20cc[5]=t20c[2]*t20c[2];


        t22aa[0]=t22a[0]*t22a[0];   t22aa[1]=t22a[1]*t22a[0];   t22aa[2]=t22a[2]*t22a[0];
                                    t22aa[3]=t22a[1]*t22a[1];   t22aa[4]=t22a[2]*t22a[1];
                                                                t22aa[5]=t22a[2]*t22a[2];

        t22bb[0]=t22b[0]*t22b[0];   t22bb[1]=t22b[1]*t22b[0];   t22bb[2]=t22b[2]*t22b[0];
                                    t22bb[3]=t22b[1]*t22b[1];   t22bb[4]=t22b[2]*t22b[1];
                                                                t22bb[5]=t22b[2]*t22b[2];

        t22cc[0]=t22c[0]*t22c[0];   t22cc[1]=t22c[1]*t22c[0];   t22cc[2]=t22c[2]*t22c[0];
                                    t22cc[3]=t22c[1]*t22c[1];   t22cc[4]=t22c[2]*t22c[1];
                                                                t22cc[5]=t22c[2]*t22c[2];

        //tensory w bazie direktora
        vect t20mol(0.0,6);
        vect t22mol(0.0,6);
        t20mol = std::sqrt(1.5)*(t20cc - Identity(3)/3.);
        t22mol = (t22aa-t22bb)/std::sqrt(2.);

	//std::cout << "sgn = " << sgn(t20es[0].e) << std::endl;

        mean_d200 = sgn(t20es[0].e)*MatrixDotProduct(t20mol,t20);
        mean_d220 = sgn(t20es[0].e)*MatrixDotProduct(t22mol,t20);
        mean_d202 = sgn(t20es[0].e)*MatrixDotProduct(t20mol,t22);
        mean_d222 = sgn(t20es[0].e)*MatrixDotProduct(t22mol,t22);

    }
    const Value & TemporalMeanEnergyPerMolecule() const {
        return energy;
    }
    const Value & TemporalMeanParity() const {
        return parity;
    }
    const Value & ParitySusceptibility() const {
        return parity_sus;
    }
    const double & MeanDelta200() const {
        return mean_d200;
    }
    const double & MeanDelta220() const {
        return mean_d220;
    }
    const double & MeanDelta202() const {
        return mean_d202;
    }
    const double & MeanDelta222() const {
        return mean_d222;
    }
    const Value & SpecificHeat() const {
        //return specific_heat;
        return fluctuation;
    }
    const Value & Fluctuation() const {
        return fluctuation;
    }


    const Value & Delta200ZByCorrelation() const {
        return d200z_from_correlation;
    }
    const Value & Delta222ZByCorrelation() const {
        return d222z_from_correlation;
    }
    const Value & Delta200XByCorrelation() const {
        return d200x_from_correlation;
    }
    const Value & Delta222XByCorrelation() const {
        return d222x_from_correlation;
    }
    const Value & Delta200YByCorrelation() const {
        return d200y_from_correlation;
    }
    const Value & Delta222YByCorrelation() const {
        return d222y_from_correlation;
    }

    const Value & Delta322ByCorrelation() const {
        return d322_from_correlation;
    }
    const Value & ParityByCorrelation() const {
        return parity_from_correlation;
    }

 
    const Value & Delta200ZByCorrelationSusceptibility() const {
        return d200z_from_correlation_sus;
    }
    const Value & Delta222ZByCorrelationSusceptibility() const {
        return d222z_from_correlation_sus;
    }
    const Value & Delta200XByCorrelationSusceptibility() const {
        return d200x_from_correlation_sus;
    }
    const Value & Delta222XByCorrelationSusceptibility() const {
        return d222x_from_correlation_sus;
    }
    const Value & Delta200YByCorrelationSusceptibility() const {
        return d200y_from_correlation_sus;
    }
    const Value & Delta222YByCorrelationSusceptibility() const {
        return d222y_from_correlation_sus;
    }

    const Value & Delta322ByCorrelationSusceptibility() const {
        return d322_from_correlation_sus;
    }
    const Value & ParityByCorrelationSusceptibility() const {
        return parity_from_correlation_sus;
    }


    const vect & Delta200ZMeanCorrelation() const {
        return mean_d200corz;
    }
    const vect & Delta222ZMeanCorrelation() const {
        return mean_d222corz;
    }
    const vect & Delta220ZMeanCorrelation() const {
        return mean_d220corz;
    }
    const vect & Delta200XMeanCorrelation() const {
        return mean_d200corx;
    }
    const vect & Delta222XMeanCorrelation() const {
        return mean_d222corx;
    }
    const vect & Delta220XMeanCorrelation() const {
        return mean_d220corx;
    }
    const vect & Delta200YMeanCorrelation() const {
        return mean_d200cory;
    }
    const vect & Delta222YMeanCorrelation() const {
        return mean_d222cory;
    }
    const vect & Delta220YMeanCorrelation() const {
        return mean_d220cory;
    }


    const vect & Delta322MeanCorrelation() const {
        return mean_d322cor;
    }
    const vect & ParityMeanCorrelation() const {
        return mean_paritycor;
    }

    const vect & MeanQxTensor() const {
        return mean_qx;
    }
    const vect & MeanQyTensor() const {
        return mean_qy;
    }
    const vect & MeanQzTensor() const {
        return mean_qz;
    }

    const double & Temperature() const {
        return temperature;
    }
    const double & Tau() const {
        return tau;
    }
    const double & Lambda() const {
        return lambda;
    }
    const double & Field() const {
        return h;
    }
};


///zapisywanie, jeżeil się doda nowe rzeczy na końcu, to nie będzie konfliktu z danymi zapisywanymi w starym formacie
template <class serializer_t>
void operator|(serializer_t & s, PRE79MeanProperties & p){
    s|p.energy;
    //s|p.specific_heat;
    s|p.fluctuation;
    s|p.d200z_from_correlation;
    s|p.d200z_from_correlation;
    s|p.d200x_from_correlation;
    s|p.d200x_from_correlation;
    s|p.d200y_from_correlation;
    s|p.d200y_from_correlation;
    s|p.d222z_from_correlation;
    s|p.d222z_from_correlation;
    s|p.d222x_from_correlation;
    s|p.d222x_from_correlation;
    s|p.d222y_from_correlation;
    s|p.d222y_from_correlation;
    s|p.d322_from_correlation;
    s|p.parity_from_correlation;
    s|p.mean_d200corz;
    s|p.mean_d220corz;
    s|p.mean_d222corz;
    s|p.mean_d200corx;
    s|p.mean_d220corx;
    s|p.mean_d222corx;
    s|p.mean_d200cory;
    s|p.mean_d220cory;
    s|p.mean_d222cory;
    s|p.mean_d322cor;
    s|p.mean_paritycor;
//    s|p.mean_d200;
//    s|p.mean_d220;
//    s|p.mean_d202;
//    s|p.mean_d222;
    s|p.mean_qx;
    s|p.mean_qy;
    s|p.mean_qz;
    s|p.temperature;
    s|p.tau;
    s|p.lambda;
    s|p.h;
    s|p.parity;
    
    s|p.parity_sus;
    s|p.d200z_from_correlation_sus;
    s|p.d200x_from_correlation_sus;
    s|p.d200y_from_correlation_sus;
    s|p.d222z_from_correlation_sus;
    s|p.d222x_from_correlation_sus;
    s|p.d222y_from_correlation_sus;
    s|p.d322_from_correlation_sus;
    s|p.parity_from_correlation_sus;
}
inline std::ostream & operator<<(std::ostream & o,const PRE79MeanProperties & p){
    o << "Temperature=" << p.Temperature() << std::endl;
    o << "Lambda=" << p.Lambda() << std::endl;
    o << "Tau=" << p.Tau() << std::endl;
    o << "Field=" << p.Field() << std::endl;
    o << "MeanEPM=" << p.TemporalMeanEnergyPerMolecule().Print() << std::endl;
    o << "MeanParity=" << p.TemporalMeanParity().Print() << std::endl;
    //o << "SpecHeat=" << p.SpecificHeat().Print() << std::endl;
    o << "Fluctuation=" << p.Fluctuation().Print() << std::endl;

    o << "Delta200ZByCorrelation=" << p.Delta200ZByCorrelation().Print() << std::endl;
    o << "Delta200XByCorrelation=" << p.Delta200XByCorrelation().Print() << std::endl;
    o << "Delta200YByCorrelation=" << p.Delta200YByCorrelation().Print() << std::endl;

    o << "Delta222ZByCorrelation=" << p.Delta222ZByCorrelation().Print() << std::endl;
    o << "Delta222XByCorrelation=" << p.Delta222XByCorrelation().Print() << std::endl;
    o << "Delta222YByCorrelation=" << p.Delta222YByCorrelation().Print() << std::endl;

    o << "TetrahedralOrder=" << p.Delta322ByCorrelation().Print() << std::endl;
    o << "ParityOrder=" << p.ParityByCorrelation().Print() << std::endl;

    o << "Delta200ZMeanCorrelation=\n" << p.Delta200ZMeanCorrelation() << std::endl;
    o << "Delta222ZMeanCorrelation=\n" << p.Delta222ZMeanCorrelation() << std::endl;
    o << "Delta220ZMeanCorrelation=\n" << p.Delta220ZMeanCorrelation() << std::endl;
    o << "Delta200XMeanCorrelation=\n" << p.Delta200XMeanCorrelation() << std::endl;
    o << "Delta222XMeanCorrelation=\n" << p.Delta222XMeanCorrelation() << std::endl;
    o << "Delta220XMeanCorrelation=\n" << p.Delta220XMeanCorrelation() << std::endl;
    o << "Delta200YMeanCorrelation=\n" << p.Delta200YMeanCorrelation() << std::endl;
    o << "Delta222YMeanCorrelation=\n" << p.Delta222YMeanCorrelation() << std::endl;
    o << "Delta220YMeanCorrelation=\n" << p.Delta220YMeanCorrelation() << std::endl;

    o << "TetrahedralMeanCorrelation=\n" << p.Delta322MeanCorrelation() << std::endl;
    o << "ParityMeanCorrelation=\n" << p.ParityMeanCorrelation() << std::endl;
    return o;
}
#endif	/* _STANDARDPROPERTIES_H */


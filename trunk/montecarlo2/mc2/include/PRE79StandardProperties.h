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
    Value   specific_heat;  ///<ciepło właściwe w funkcji czasu (w praktyce obliczanie całej historii jest bardzo kosztowne, obliczamy dopiero na końcu)
    Value   fluctuation;

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
public:
    ///obliczanie ciepła właściwego metodą bootstrapu
    void CalculateSpecificHeat(const double & T){
        /*
        vect fluct(acc_idx+1);
        for(int i=0;i<(acc_idx+1);i++){
            double E = BootstrapMean(energy,0,acc_idx+1,1);
            double E2 = BootstrapMean(vect(energy*energy),0,acc_idx+1,1);
            fluct[i] = (E2-E*E)*lat->GetN()/T/T;
        }
        //specific_heat=BootstrapMean(fluct);
        specific_heat=Mean(fluct);
         */

        
        //std::cout << "TEST SPH from Boostrap/Boostrap: " << specific_heat << std::endl;
        //std::cout << "TEST SPH from Boostrap/Mean: " << Mean(fluct) << std::endl;

        vect fluct2(acc_idx+1);
        vect fluct3(acc_idx+1);
        for(int i=0;i<(acc_idx+1);i++){
            vect e(acc_idx+1);
            vect e2(acc_idx+1);
            for(int j=0;j<(acc_idx+1);j++){
                int t = (acc_idx)*random01();
                e[j]=energy[t];
                e2[j]=energy[t]*energy[t];
            }
            fluct2[i]=(double(Mean(e2))-double(Mean(e))*double(Mean(e)))*lat->GetN()/T/T;
            fluct3[i]=(double(Mean(e2))-double(Mean(e))*double(Mean(e)))*lat->GetN();
        }
        //std::cout << "TEST SPH from old BS: " << BootstrapMean(fluct2) << std::endl;
        specific_heat=BootstrapMean(fluct2);
        fluctuation=BootstrapMean(fluct3);

        //std::cout << "TEST SPH from Fluctuation by bootstrap: " << (std::pow(BootstrapMean(energy).Error(),2.0)*lat->GetN()/T/T) << std::endl;
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
    paritycor(l,_ncycles)
    {
        acc_idx=-1;
        index=0;
        energy.resize(ncycles,0.0);
        for(int i=0;i<ncycles;i++){
            MeanQxTensor[i].resize(6,0.0);
            MeanQyTensor[i].resize(6,0.0);
            MeanQzTensor[i].resize(6,0.0);
        }
        //specific_heat.resize(ncycles);
    }


    PRE79StandardProperties(const PRE79StandardProperties & p){
        //tak to musi być zrobione, por. Standard C++ Library Reference pp. 327
        energy.resize(p.energy.size());
        energy = p.energy;
        lat=p.lat;
        readonly=p.readonly;
        index=p.index;
        ncycles=p.ncycles;
        specific_heat=p.specific_heat;
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
        energy.resize(p.energy.size());
        energy = p.energy;
        lat=p.lat;
        readonly=p.readonly;
        index=p.index;
        ncycles=p.ncycles;
        specific_heat=p.specific_heat;
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
        CalculateMeanTensors();

        if((acc_idx+1)>=ncycles)
            CalculateSpecificHeat(H->GetTemperature());
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
        for(int i=0;i<oldsize;i++)
            tmpenergy[i]=energy[i];

        for(int i=oldsize;i<ncycles;i++){
            tmpenergy[i]=p.energy[i-oldsize];
        }
        energy.resize(ncycles,0.0);
        energy=tmpenergy;
    }

    Value TemporalMeanEnergyPerMolecule() const {
        return BootstrapMean(energy,0,acc_idx+1);
    }

    //korelacje dla różnych osi
    Value Delta200ZByCorrelation() const {
        return sqrt(d200corz.Limit());
    }
    Value Delta222ZByCorrelation() const {
        return sqrt(d222corz.Limit());
    }
    Value Delta200XByCorrelation() const {
        return sqrt(d200corx.Limit());
    }
    Value Delta222XByCorrelation() const {
        return sqrt(d222corx.Limit());
    }
    Value Delta200YByCorrelation() const {
        return sqrt(d200cory.Limit());
    }
    Value Delta222YByCorrelation() const {
        return sqrt(d222cory.Limit());
    }


    Value Delta322ByCorrelation() const {
        return sqrt(d322cor.Limit());
    }
    Value ParityByCorrelation() const {
        return sqrt(paritycor.Limit());
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
        return specific_heat;
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
    s|prop.specific_heat;
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
    s|prop.d322cor;
    s|prop.paritycor;
    s|prop.MeanQxTensor;
    s|prop.MeanQyTensor;
    s|prop.MeanQzTensor;
    s|prop.energy;
}

/**
 * używane do przechowania końcowych, uśrednionych parametrów układu
 */
class PRE79MeanProperties {
    template<class serializer_t>
    friend void operator|(serializer_t & s, PRE79MeanProperties & prop);
    Value energy;
    Value specific_heat;
    Value fluctuation;
    Value d200z_from_correlation;
    Value d222z_from_correlation;
    Value d200x_from_correlation;
    Value d222x_from_correlation;
    Value d200y_from_correlation;
    Value d222y_from_correlation;
    Value d322_from_correlation;
    Value parity_from_correlation;

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
        specific_heat = prop.SpecificHeat();
        fluctuation = prop.Fluctuation();

        d200z_from_correlation = prop.Delta200ZByCorrelation();
        d222z_from_correlation = prop.Delta222ZByCorrelation();
        d200x_from_correlation = prop.Delta200XByCorrelation();
        d222x_from_correlation = prop.Delta222XByCorrelation();
        d200y_from_correlation = prop.Delta200YByCorrelation();
        d222y_from_correlation = prop.Delta222YByCorrelation();


        d322_from_correlation = prop.Delta322ByCorrelation();
        parity_from_correlation = prop.ParityByCorrelation();

        mean_d200corz = prop.Delta200ZMeanCorrelation();
        mean_d220corz = prop.Delta220ZMeanCorrelation();
        mean_d222corz = prop.Delta222ZMeanCorrelation();
        mean_d200corx = prop.Delta200XMeanCorrelation();
        mean_d220corx = prop.Delta220XMeanCorrelation();
        mean_d222corx = prop.Delta222XMeanCorrelation();
        mean_d200cory = prop.Delta200YMeanCorrelation();
        mean_d220cory = prop.Delta220YMeanCorrelation();
        mean_d222cory = prop.Delta222YMeanCorrelation();

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
    const Value & TemporalMeanEnergyPerMolecule() const {
        return energy;
    }
    const Value & SpecificHeat() const {
        return specific_heat;
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

template <class serializer_t>
void operator|(serializer_t & s, PRE79MeanProperties & p){
    s|p.energy;
    s|p.specific_heat;
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
    s|p.mean_qx;
    s|p.mean_qy;
    s|p.mean_qz;
    s|p.temperature;
    s|p.tau;
    s|p.lambda;
    s|p.h;
}
inline std::ostream & operator<<(std::ostream & o,const PRE79MeanProperties & p){
    o << "Temperature=" << p.Temperature() << std::endl;
    o << "Lambda=" << p.Lambda() << std::endl;
    o << "Tau=" << p.Tau() << std::endl;
    o << "Field=" << p.Field() << std::endl;
    o << "MeanEPM=" << p.TemporalMeanEnergyPerMolecule().Print() << std::endl;
    o << "SpecHeat=" << p.SpecificHeat().Print() << std::endl;
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
}
#endif	/* _STANDARDPROPERTIES_H */


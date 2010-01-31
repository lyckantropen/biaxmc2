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
    Value    specific_heat;  ///<ciepło właściwe w funkcji czasu (w praktyce obliczanie całej historii jest bardzo kosztowne, obliczamy dopiero na końcu)

    //funkcje korelacji dla różnych czasów
    Delta200Correlation d200cor;
    Delta222Correlation d222cor;
    Delta220Correlation d220cor;
    Delta322Correlation d322cor;
    ParityCorrelation paritycor;

    ///obliczanie średniej energii (po siatce)
    void CalculateMeanEnergy(){
        double E=0.0;
        for(int site=0;site<lat->GetN();site++){
            E+=lat->Particles[site].GetEnergy();
        }
        E/=double(lat->GetN());
        energy[acc_idx]=E;
    }
    ///obliczanie ciepła właściwego metodą bootstrapu
    void CalculateSpecificHeat(const double & T){
        vect fluct(acc_idx+1);
        for(int i=0;i<(acc_idx+1);i++){
            double E = BootstrapMean(energy,0,acc_idx+1);
            double E2 = BootstrapMean(vect(energy*energy),0,acc_idx+1);
            fluct[i] = (E2-E*E)*lat->GetN()/T/T;
        }
        specific_heat=BootstrapMean(fluct);
    }
public:
    ///konstruktor serializacyjny
    PRE79StandardProperties():readonly(true),index(0),acc_idx(0) {}
    ///konstruktor tradycyjny
    PRE79StandardProperties(Lattice * l, int _ncycles):
    lat(l),readonly(false),ncycles(_ncycles),
    d200cor(l,_ncycles),
    d222cor(l,_ncycles),
    d220cor(l,_ncycles),
    d322cor(l,_ncycles),
    paritycor(l,_ncycles)
    {
        acc_idx=-1;
        index=0;
        energy.resize(ncycles,0.0);
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
        d200cor=p.d200cor;
        d220cor=p.d220cor;
        d222cor=p.d222cor;
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
        d200cor=p.d200cor;
        d220cor=p.d220cor;
        d222cor=p.d222cor;
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

        d200cor.Update();
        d220cor.Update();
        d222cor.Update();
        d322cor.Update();
        paritycor.Update();

        if((acc_idx+1)>=ncycles)
            CalculateSpecificHeat(H->GetTemperature());
        CalculateMeanEnergy();
    }
    Value TemporalMeanEnergyPerMolecule() const {
        return BootstrapMean(energy,0,acc_idx+1);
    }
    Value UniaxialOrderByCorrelation() const {
        return sqrt(d200cor.Limit());
    }
    Value BiaxialOrderByCorrelation() const {
        return sqrt(d222cor.Limit());
    }
    Value TetrahedralOrderByCorrelation() const {
        return sqrt(d322cor.Limit());
    }
    Value ParityOrderByCorrelation() const {
        return sqrt(paritycor.Limit());
    }
    vect UniaxialMeanCorrelation() const {
        return d200cor.Mean();
    }
    vect BiaxialMeanCorrelation() const {
        return d222cor.Mean();
    }
    vect TetrahedralMeanCorrelation() const {
        return d322cor.Mean();
    }
    vect ParityMeanCorrelation() const {
        return paritycor.Mean();
    }
    vect Delta220MeanCorrelation() const {
        return d220cor.Mean();
    }
    const vect & GetUniaxialCorrelation() const {
        return d200cor.CurrentCorrelation();
    }
    const vect & GetBiaxialCorrelation() const {
        return d222cor.CurrentCorrelation();
    }
    const vect & GetTetrahedralCorrelation() const {
        return d322cor.CurrentCorrelation();
    }
    const vect & GetParityCorrelation() const {
        return paritycor.CurrentCorrelation();
    }
    const vect & GetDelta220Correlation() const {
        return d220cor.CurrentCorrelation();
    }
    const vect & EnergyEvolution() const {
        return energy;
    }
    const Delta200Correlation & UniaxialCorrelationEvolution() const {
        return d200cor;
    }
    const Delta222Correlation & BiaxialCorrelationEvolution() const {
        return d222cor;
    }
    const Delta322Correlation & TetrahedralCorrelationEvolution() const {
        return d322cor;
    }
    const Delta220Correlation & Delta220CorrelationEvolution() const {
        return d220cor;
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
    const int & GetIndex() const {
        return index;
    }
    const int & GetAccIdx() const {
        return acc_idx;
    }
    const int & GetNCycles() const {
        return ncycles;
    }

};

template<class serializer_t>
void operator|(serializer_t & s, PRE79StandardProperties & prop){
    s|prop.index;
    s|prop.acc_idx;
    s|prop.ncycles;
    s|prop.specific_heat;
    s|prop.d200cor;
    s|prop.d220cor;
    s|prop.d222cor;
    s|prop.d322cor;
    s|prop.paritycor;
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
    Value uniaxial_from_correlation;
    Value biaxial_from_correlation;
    Value tetrahedral_from_correlation;
    Value parity_from_correlation;

    vect mean_d200cor;
    vect mean_d220cor;
    vect mean_d222cor;
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
        energy = prop.TemporalMeanEnergyPerMolecule();
        specific_heat = prop.SpecificHeat();
        uniaxial_from_correlation = prop.UniaxialOrderByCorrelation();
        biaxial_from_correlation = prop.BiaxialOrderByCorrelation();
        tetrahedral_from_correlation = prop.TetrahedralOrderByCorrelation();
        parity_from_correlation = prop.ParityOrderByCorrelation();

        mean_d200cor = prop.UniaxialMeanCorrelation();
        mean_d220cor = prop.Delta220MeanCorrelation();
        mean_d222cor = prop.BiaxialMeanCorrelation();
        mean_d322cor = prop.TetrahedralMeanCorrelation();
        mean_paritycor = prop.ParityMeanCorrelation();

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
    const Value & UniaxialOrderByCorrelation() const {
        return uniaxial_from_correlation;
    }
    const Value & BiaxialOrderByCorrelation() const {
        return biaxial_from_correlation;
    }
    const Value & TetrahedralOrderByCorrelation() const {
        return tetrahedral_from_correlation;
    }
    const Value & ParityOrderByCorrelation() const {
        return parity_from_correlation;
    }

    const vect & UniaxialMeanCorrelation() const {
        return mean_d200cor;
    }
    const vect & BiaxialMeanCorrelation() const {
        return mean_d222cor;
    }
    const vect & Delta220MeanCorrelation() const {
        return mean_d220cor;
    }
    const vect & TetrahedralMeanCorrelation() const {
        return mean_d322cor;
    }
    const vect & ParityMeanCorrelation() const {
        return mean_paritycor;
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
    s|p.uniaxial_from_correlation;
    s|p.biaxial_from_correlation;
    s|p.tetrahedral_from_correlation;
    s|p.parity_from_correlation;
    s|p.mean_d200cor;
    s|p.mean_d220cor;
    s|p.mean_d322cor;
    s|p.mean_paritycor;
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
    o << "UniaxialOrder=" << p.UniaxialOrderByCorrelation().Print() << std::endl;
    o << "BiaxialOrder=" << p.BiaxialOrderByCorrelation().Print() << std::endl;
    o << "TetrahedralOrder=" << p.TetrahedralOrderByCorrelation().Print() << std::endl;
    o << "ParityOrder=" << p.ParityOrderByCorrelation().Print() << std::endl;
    o << "UniaxialMeanCorrelation=\n" << p.UniaxialMeanCorrelation() << std::endl;
    o << "BiaxialMeanCorrelation=\n" << p.BiaxialMeanCorrelation() << std::endl;
    o << "Delta220MeanCorrelation=\n" << p.Delta220MeanCorrelation() << std::endl;
    o << "TetrahedralMeanCorrelation=\n" << p.TetrahedralMeanCorrelation() << std::endl;
    o << "ParityMeanCorrelation=\n" << p.ParityMeanCorrelation() << std::endl;
}
#endif	/* _STANDARDPROPERTIES_H */


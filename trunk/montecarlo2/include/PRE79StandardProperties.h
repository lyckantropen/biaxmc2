/*
 * File:   StandardProperties.h
 * Author: karol
 *
 * Created on 18 listopad 2009, 11:27
 */

#ifndef _STANDARDPROPERTIES_H
#define _STANDARDPROPERTIES_H

#include "SpatialCorrelationEvolution.h"
#include "PRE79StandardHamiltonian.h"
#include "AutoCorrelationTimeCalculator.h"
#include "Lattice.h"
#include "Statistical.h"
#include "Contractions.h"
#include "PRE79SpatialCorrelations.h"

#include "serializer.h"
#include "evsort.h"


/**
 * @brief Measure physical quantities and store their history
 *
 * This class:
 * - collects the properties of the system during measurement
 * - stores the entire history of properties in vectors
 * - calculates some of the quantities (e.g. specific heat)
 *
 */
class PRE79StandardProperties
{
    ///Serializer operator
    template<class serializer_t>
    friend void operator|(serializer_t & s, PRE79StandardProperties & prop);

    const Lattice *   lat;        ///<pointer to the lattice we are operating on
    bool    readonly;       ///<readonly flag
    int index;              ///<unused
    int acc_idx;            ///<number of the present measurement
    int ncycles;            ///<total number of measurements

    vect    energy;         ///<mean energy per molecule as a function of time
    vect    parity;         ///<mean parity per molecule as a function
    vect    T20T20z;         ///<T20z*T20z as a function of time
    vect    T22T22z;         ///<T22z*T22z as a function of time
    vect    T20T20x;         ///<T20x*T20x as a function of time
    vect    T22T22x;         ///<T22x*T22x as a function of time
    vect    T20T20y;         ///<T20y*T20y as a function of time
    vect    T22T22y;         ///<T22y*T22y as a function of time
    vect    T32T32;          ///<T32*T32 as a function of time

    vect        autocorrelation_time;   ///<autocorrelation time as a function of time

    Value   fluctuation;    ///<fluctuation of energy, calculated at the end


    //ewolucja średnich tensorów xx,yy,zz
    std::vector<vect>  MeanQxTensor;            ///<history of the mean tensor <X(x)X>
    std::vector<vect>  MeanQyTensor;            ///<history of the mean tensor <Y(x)Y>
    std::vector<vect>  MeanQzTensor;            ///<history of the mean tensor <Z(x)Z>

    //ewolucja funkcji korelacji dla różnych czasów
    Delta200CorrelationZ d200corz;      ///<correlation function history for T20z(0)*T20z(r)
    Delta222CorrelationZ d222corz;      ///<correlation function history for T22z(0)*T22z(r)
    Delta220CorrelationZ d220corz;      ///<correlation function history for T22z(0)*T20z(r)

    Delta200CorrelationX d200corx;      ///<correlation function history for T20x(0)*T20x(r)
    Delta222CorrelationX d222corx;      ///<correlation function history for T22x(0)*T22x(r)
    Delta220CorrelationX d220corx;      ///<correlation function history for T22x(0)*T20x(r)

    Delta200CorrelationY d200cory;      ///<correlation function history for T20y(0)*T20y(r)
    Delta222CorrelationY d222cory;      ///<correlation function history for T22y(0)*T22y(r)
    Delta220CorrelationY d220cory;      ///<correlation function history for T22y(0)*T20y(r)

    Delta322Correlation d322cor; ///<correlation function history for T32(0)*T32(r)
    ParityCorrelation paritycor; ///<correlation function history for p(0)*p(r)

    ///bookkeeping of inverse temperature indices for parallel tempering
    vect b_idx;
    ///bookkeeping of inverse temperature values for parallel tempering
    vect b_val;
    ///current beta value for parallel tempering
    double cur_b_val;
    ///current index value for parallel tempering
    double cur_b_idx;

    /**
     * Calculate the instanteneous lattice mean energy and update the energy history vector.
     */
    void CalculateMeanEnergy();
    /**
     * Calculate the instanteneous lattice mean parity and update the parity history vector
     */
    void CalculateMeanParity();
public:
    /**
     * Calculate the fluctuation of energy per molecule. To obtain specific heat,
     * divide by T^2
     */
    void CalculateSpecificHeat();
private:
    ///obliczanie średnich tensorów xx,yy i zz
    void CalculateMeanTensors();

public:
    ///konstruktor serializacyjny
    PRE79StandardProperties();
    ///konstruktor tradycyjny
    PRE79StandardProperties(const Lattice * l, int _ncycles);

    ///obliczenie kontrakcji
    void Update(int idx, const StandardL2Hamiltonian * H,
                const AutoCorrelationTimeCalculator * ac = NULL);

    ///dołożenie danych z innych właściwości, koniecznie po zakończeniu obliczeń
    void Append(const PRE79StandardProperties * _p);

    Value TemporalMeanEnergyPerMolecule() const
    {
        //return BootstrapMean(energy,0,acc_idx+1);
        return Mean(energy, 0, acc_idx + 1);
    }
    Value TemporalMeanParity() const
    {
        return Mean(parity, 0, acc_idx + 1);
    }
    Value ParitySusceptibility() const
    {
        return CalculateFluctuation(parity, acc_idx);
    }
    //korelacje dla różnych osi
    Value Delta200ZByCorrelation() const
    {
        return sqrt(d200corz.Limit());
    }
    Value Delta200ZByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d200corz.LimitHistory())));
    }
    Value Delta222ZByCorrelation() const
    {
        return sqrt(d222corz.Limit());
    }
    Value Delta222ZByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d222corz.LimitHistory())));
    }
    Value Delta200XByCorrelation() const
    {
        return sqrt(d200corx.Limit());
    }
    Value Delta200XByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d200corx.LimitHistory())));
    }
    Value Delta222XByCorrelation() const
    {
        return sqrt(d222corx.Limit());
    }
    Value Delta222XByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d222corx.LimitHistory())));
    }
    Value Delta200YByCorrelation() const
    {
        return sqrt(d200cory.Limit());
    }
    Value Delta200YByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d200cory.LimitHistory())));
    }
    Value Delta222YByCorrelation() const
    {
        return sqrt(d222cory.Limit());
    }
    Value Delta222YByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d222cory.LimitHistory())));
    }

    Value Delta322ByCorrelation() const
    {
        return sqrt(d322cor.Limit());
    }
    Value Delta322ByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(d322cor.LimitHistory())));
    }
    Value ParityByCorrelation() const
    {
        return sqrt(paritycor.Limit());
    }
    Value ParityByCorrelationSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(std::abs(paritycor.LimitHistory())));
    }

    Value MeanT20T20Z() const
    {
        return Mean(T20T20z, 0, acc_idx + 1);
    }
    Value MeanT20T20ZSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T20T20z));
    }

    Value MeanT22T22Z() const
    {
        return Mean(T22T22z, 0, acc_idx + 1);
    }
    Value MeanT22T22ZSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T22T22z));
    }
    Value MeanT20T20X() const
    {
        return Mean(T20T20x, 0, acc_idx + 1);
    }
    Value MeanT20T20XSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T20T20x));
    }

    Value MeanT22T22X() const
    {
        return Mean(T22T22x, 0, acc_idx + 1);
    }
    Value MeanT22T22XSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T22T22x));
    }
    Value MeanT20T20Y() const
    {
        return Mean(T20T20y, 0, acc_idx + 1);
    }
    Value MeanT20T20YSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T20T20y));
    }

    Value MeanT22T22Y() const
    {
        return Mean(T22T22y, 0, acc_idx + 1);
    }
    Value MeanT22T22YSusceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T22T22y));
    }

    Value MeanT32T32() const
    {
        return Mean(T32T32, 0, acc_idx + 1);
    }
    Value MeanT32T32Susceptibility() const
    {
        return CalculateFluctuation(std::sqrt(T32T32));
    }

    Value MeanAutocorrelationTime() const
    {
        return Mean(autocorrelation_time, 0, acc_idx + 1);
    }


    //średnie funkcje korelacji dla poszczególnych osi
    vect Delta200ZMeanCorrelation() const
    {
//        std::cout << "D200Z\n";
        return d200corz.Mean(acc_idx + 1);
    }
    vect Delta222ZMeanCorrelation() const
    {
//        std::cout << "D222Z\n";
        return d222corz.Mean(acc_idx + 1);
    }
    vect Delta220ZMeanCorrelation() const
    {
//        std::cout << "D220Z\n";
        return d220corz.Mean(acc_idx + 1);
    }
    vect Delta200YMeanCorrelation() const
    {
//        std::cout << "D200Y\n";
        return d200cory.Mean(acc_idx + 1);
    }
    vect Delta222YMeanCorrelation() const
    {
//        std::cout << "D222Y\n";
        return d222cory.Mean(acc_idx + 1);
    }
    vect Delta220YMeanCorrelation() const
    {
//        std::cout << "D220Y\n";
        return d220cory.Mean(acc_idx + 1);
    }
    vect Delta200XMeanCorrelation() const
    {
//        std::cout << "D200X\n";
        return d200corx.Mean(acc_idx + 1);
    }
    vect Delta222XMeanCorrelation() const
    {
//        std::cout << "D222X\n";
        return d222corx.Mean(acc_idx + 1);
    }
    vect Delta220XMeanCorrelation() const
    {
//        std::cout << "D220X\n";
        return d220corx.Mean(acc_idx + 1);
    }



    vect Delta322MeanCorrelation() const
    {
//        std::cout << "D322\n";
        return d322cor.Mean(acc_idx + 1);
    }
    vect ParityMeanCorrelation() const
    {
//        std::cout << "Parity\n";
        return paritycor.Mean(acc_idx + 1);
    }


    const vect & GetDelta200ZCorrelation() const
    {
        return d200corz.CurrentCorrelation();
    }
    const vect & GetDelta222ZCorrelation() const
    {
        return d222corz.CurrentCorrelation();
    }
    const vect & GetDelta220ZCorrelation() const
    {
        return d220corz.CurrentCorrelation();
    }
    const vect & GetDelta200XCorrelation() const
    {
        return d200corx.CurrentCorrelation();
    }
    const vect & GetDelta222XCorrelation() const
    {
        return d222corx.CurrentCorrelation();
    }
    const vect & GetDelta220XCorrelation() const
    {
        return d220corx.CurrentCorrelation();
    }
    const vect & GetDelta200YCorrelation() const
    {
        return d200cory.CurrentCorrelation();
    }
    const vect & GetDelta222YCorrelation() const
    {
        return d222cory.CurrentCorrelation();
    }
    const vect & GetDelta220YCorrelation() const
    {
        return d220cory.CurrentCorrelation();
    }


    const vect & GetDelta322Correlation() const
    {
        return d322cor.CurrentCorrelation();
    }
    const vect & GetParityCorrelation() const
    {
        return paritycor.CurrentCorrelation();
    }

    vect GetMeanQxTensor() const
    {
//        std::cout << "Qx\n";
        return MeanVector(MeanQxTensor, 0, acc_idx + 1);
    }
    vect GetMeanQyTensor() const
    {
//        std::cout << "Qy\n";
        return MeanVector(MeanQyTensor, 0, acc_idx + 1);
    }
    vect GetMeanQzTensor() const
    {
//        std::cout << "Qz\n";
        return MeanVector(MeanQzTensor, 0, acc_idx + 1);
    }

    const vect & EnergyEvolution() const
    {
        return energy;
    }
    const vect & ParityEvolution() const
    {
        return parity;
    }

    const Delta200CorrelationZ & Delta200ZCorrelationEvolution() const
    {
        return d200corz;
    }
    const Delta222CorrelationZ & Delta222ZCorrelationEvolution() const
    {
        return d222corz;
    }
    const Delta220CorrelationZ & Delta220ZCorrelationEvolution() const
    {
        return d220corz;
    }
    const Delta200CorrelationX & Delta200XCorrelationEvolution() const
    {
        return d200corx;
    }
    const Delta222CorrelationX & Delta222XCorrelationEvolution() const
    {
        return d222corx;
    }
    const Delta220CorrelationX & Delta220XCorrelationEvolution() const
    {
        return d220corx;
    }
    const Delta200CorrelationY & Delta200YCorrelationEvolution() const
    {
        return d200cory;
    }
    const Delta222CorrelationY & Delta222YCorrelationEvolution() const
    {
        return d222cory;
    }
    const Delta220CorrelationY & Delta220YCorrelationEvolution() const
    {
        return d220cory;
    }

    const Delta322Correlation & Delta322CorrelationEvolution() const
    {
        return d322cor;
    }

    const ParityCorrelation & ParityCorrelationEvolution() const
    {
        return paritycor;
    }

    //Accessors
    ///ciepło właściwe
    const Value & SpecificHeat() const
    {
        //CalculateSpecificHeat(H->GetTemperature());
        //return specific_heat;
        return fluctuation;
    }
    const vect & GetAutocorrelationTimeHistory() const
    {
        return autocorrelation_time;
    }

    const Value & Fluctuation() const
    {
        return fluctuation;
    }
    const int & GetIndex() const
    {
        return index;
    }
    const int & GetAccIdx() const
    {
        return acc_idx;
    }
    const int & GetNCycles() const
    {
        return ncycles;
    }
    const vect & GetBetaValues() const
    {
        return b_val;
    }
    const vect & GetBetaIndices() const
    {
        return b_idx;
    }
    int GetMaxCorrLen() const
    {
        //return lat->GetL()/2+1;
        return paritycor.GetMax() + 1;
    }
    void SetLattice(Lattice * l)
    {
        lat = l;
    }
    const Lattice * GetLattice() const
    {
        return lat;
    }
    void SetCurrentBetaValue(const double & bv)
    {
        cur_b_val = bv;
        b_val[acc_idx] = bv;
    }
    void SetCurrentBetaIndex(const double & bi)
    {
        cur_b_idx = bi;
        b_idx[acc_idx] = bi;
    }
    const double & GetCurrentBetaValue() const
    {
        return cur_b_val;
    }
    const double & GetCurrentBetaIndex() const
    {
        return cur_b_idx;
    }
    std::valarray<bool>      GetReplicaMask(const double & replica) const;
    static PRE79StandardProperties * FromReplicaMask(const PRE79StandardProperties& p, const double & replica);
};

template<class serializer_t>
void operator|(serializer_t & s, PRE79StandardProperties & prop)
{
    s | prop.index;
    s | prop.acc_idx;
    s | prop.ncycles;
    //s|prop.specific_heat;
    s | prop.fluctuation;
    s | prop.d200corz;
    s | prop.d220corz;
    s | prop.d222corz;
    s | prop.d200corx;
    s | prop.d220corx;
    s | prop.d222corx;
    s | prop.d200cory;
    s | prop.d220cory;
    s | prop.d222cory;
//    s|prop.D200;
//    s|prop.D220;
//    s|prop.D202;
//    s|prop.D222;
    s | prop.d322cor;
    s | prop.paritycor;
    s | prop.MeanQxTensor;
    s | prop.MeanQyTensor;
    s | prop.MeanQzTensor;
    s | prop.energy;
    s | prop.parity;

    s | prop.T20T20z;
    s | prop.T22T22z;

    s | prop.T20T20x;
    s | prop.T22T22x;

    s | prop.T20T20y;
    s | prop.T22T22y;
    s | prop.T32T32;

    s | prop.autocorrelation_time;

    s | prop.cur_b_val;
    s | prop.cur_b_idx;
    s | prop.b_val;
    s | prop.b_idx;
}

/**
 * final, averaged statistical properties of the system
 */
class PRE79MeanProperties
{
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
    Value T20T20z;
    Value T22T22z;
    Value T20T20z_sus;
    Value T22T22z_sus;
    Value T20T20x;
    Value T22T22x;
    Value T20T20x_sus;
    Value T22T22x_sus;
    Value T20T20y;
    Value T22T22y;
    Value T20T20y_sus;
    Value T22T22y_sus;
    Value T32T32;
    Value T32T32_sus;
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
    double kappa;

    Value mean_autocorrelation_time;



public:
    ///konstruktor serializacyjny
    PRE79MeanProperties();
    ///konstruktor tradycyjny
    PRE79MeanProperties(PRE79StandardProperties * _prop, StandardL2Hamiltonian * _H);
//    PRE79MeanProperties(const PRE79MeanProperties & s);
//    const PRE79MeanProperties & operator=(const PRE79MeanProperties & s);
    //do wykonania po wczytaniu
    void CalculateMeanTensors();
    const Value & TemporalMeanEnergyPerMolecule() const
    {
        return energy;
    }
    const Value & TemporalMeanParity() const
    {
        return parity;
    }
    const Value & ParitySusceptibility() const
    {
        return parity_sus;
    }
    const double & MeanDelta200() const
    {
        return mean_d200;
    }
    const double & MeanDelta220() const
    {
        return mean_d220;
    }
    const double & MeanDelta202() const
    {
        return mean_d202;
    }
    const double & MeanDelta222() const
    {
        return mean_d222;
    }
    const Value & SpecificHeat() const
    {
        //return specific_heat;
        return fluctuation;
    }
    const Value & Fluctuation() const
    {
        return fluctuation;
    }


    const Value & Delta200ZByCorrelation() const
    {
        return d200z_from_correlation;
    }
    const Value & Delta222ZByCorrelation() const
    {
        return d222z_from_correlation;
    }
    const Value & Delta200XByCorrelation() const
    {
        return d200x_from_correlation;
    }
    const Value & Delta222XByCorrelation() const
    {
        return d222x_from_correlation;
    }
    const Value & Delta200YByCorrelation() const
    {
        return d200y_from_correlation;
    }
    const Value & Delta222YByCorrelation() const
    {
        return d222y_from_correlation;
    }

    const Value & Delta322ByCorrelation() const
    {
        return d322_from_correlation;
    }
    const Value & ParityByCorrelation() const
    {
        return parity_from_correlation;
    }


    const Value & Delta200ZByCorrelationSusceptibility() const
    {
        return d200z_from_correlation_sus;
    }
    const Value & Delta222ZByCorrelationSusceptibility() const
    {
        return d222z_from_correlation_sus;
    }
    const Value & Delta200XByCorrelationSusceptibility() const
    {
        return d200x_from_correlation_sus;
    }
    const Value & Delta222XByCorrelationSusceptibility() const
    {
        return d222x_from_correlation_sus;
    }
    const Value & Delta200YByCorrelationSusceptibility() const
    {
        return d200y_from_correlation_sus;
    }
    const Value & Delta222YByCorrelationSusceptibility() const
    {
        return d222y_from_correlation_sus;
    }

    const Value & Delta322ByCorrelationSusceptibility() const
    {
        return d322_from_correlation_sus;
    }
    const Value & ParityByCorrelationSusceptibility() const
    {
        return parity_from_correlation_sus;
    }
    Value MeanT20T20Z() const
    {
        return T20T20z;
    }
    Value MeanT20T20ZSusceptibility() const
    {
        return T20T20z_sus;
    }

    Value MeanT22T22Z() const
    {
        return T22T22z;
    }
    Value MeanT22T22ZSusceptibility() const
    {
        return T22T22z_sus;
    }
    Value MeanT20T20X() const
    {
        return T20T20x;
    }
    Value MeanT20T20XSusceptibility() const
    {
        return T20T20x_sus;
    }

    Value MeanT22T22X() const
    {
        return T22T22x;
    }
    Value MeanT22T22XSusceptibility() const
    {
        return T22T22x_sus;
    }
    Value MeanT20T20Y() const
    {
        return T20T20y;
    }
    Value MeanT20T20YSusceptibility() const
    {
        return T20T20y_sus;
    }

    Value MeanT22T22Y() const
    {
        return T22T22y;
    }
    Value MeanT22T22YSusceptibility() const
    {
        return T22T22y_sus;
    }

    Value MeanT32T32() const
    {
        return T32T32;
    }
    Value MeanT32T32Susceptibility() const
    {
        return T32T32_sus;
    }

    const vect & Delta200ZMeanCorrelation() const
    {
        return mean_d200corz;
    }
    const vect & Delta222ZMeanCorrelation() const
    {
        return mean_d222corz;
    }
    const vect & Delta220ZMeanCorrelation() const
    {
        return mean_d220corz;
    }
    const vect & Delta200XMeanCorrelation() const
    {
        return mean_d200corx;
    }
    const vect & Delta222XMeanCorrelation() const
    {
        return mean_d222corx;
    }
    const vect & Delta220XMeanCorrelation() const
    {
        return mean_d220corx;
    }
    const vect & Delta200YMeanCorrelation() const
    {
        return mean_d200cory;
    }
    const vect & Delta222YMeanCorrelation() const
    {
        return mean_d222cory;
    }
    const vect & Delta220YMeanCorrelation() const
    {
        return mean_d220cory;
    }


    const vect & Delta322MeanCorrelation() const
    {
        return mean_d322cor;
    }
    const vect & ParityMeanCorrelation() const
    {
        return mean_paritycor;
    }

    const vect & MeanQxTensor() const
    {
        return mean_qx;
    }
    const vect & MeanQyTensor() const
    {
        return mean_qy;
    }
    const vect & MeanQzTensor() const
    {
        return mean_qz;
    }

    const double & Temperature() const
    {
        return temperature;
    }
    const double & Tau() const
    {
        return tau;
    }
    const double & Lambda() const
    {
        return lambda;
    }
    const double & Field() const
    {
        return h;
    }
    const double & Kappa() const
    {
        return kappa;
    }

    const Value & MeanAutocorrelationTime() const
    {
        return mean_autocorrelation_time;
    }

};


///zapisywanie, jeżeil się doda nowe rzeczy na końcu, to nie będzie konfliktu z danymi zapisywanymi w starym formacie
template <class serializer_t>
void operator|(serializer_t & s, PRE79MeanProperties & p)
{
    s | p.energy;
    //s|p.specific_heat;
    s | p.fluctuation;
    s | p.d200z_from_correlation;
    s | p.d200z_from_correlation;
    s | p.d200x_from_correlation;
    s | p.d200x_from_correlation;
    s | p.d200y_from_correlation;
    s | p.d200y_from_correlation;
    s | p.d222z_from_correlation;
    s | p.d222z_from_correlation;
    s | p.d222x_from_correlation;
    s | p.d222x_from_correlation;
    s | p.d222y_from_correlation;
    s | p.d222y_from_correlation;
    s | p.d322_from_correlation;
    s | p.parity_from_correlation;
    s | p.mean_d200corz;
    s | p.mean_d220corz;
    s | p.mean_d222corz;
    s | p.mean_d200corx;
    s | p.mean_d220corx;
    s | p.mean_d222corx;
    s | p.mean_d200cory;
    s | p.mean_d220cory;
    s | p.mean_d222cory;
    s | p.mean_d322cor;
    s | p.mean_paritycor;
//    s|p.mean_d200;
//    s|p.mean_d220;
//    s|p.mean_d202;
//    s|p.mean_d222;
    s | p.mean_qx;
    s | p.mean_qy;
    s | p.mean_qz;
    s | p.temperature;
    s | p.tau;
    s | p.lambda;
    s | p.h;
    s | p.parity;

    s | p.parity_sus;
    s | p.d200z_from_correlation_sus;
    s | p.d200x_from_correlation_sus;
    s | p.d200y_from_correlation_sus;
    s | p.d222z_from_correlation_sus;
    s | p.d222x_from_correlation_sus;
    s | p.d222y_from_correlation_sus;
    s | p.d322_from_correlation_sus;
    s | p.parity_from_correlation_sus;

    s | p.T20T20z;
    s | p.T20T20z_sus;
    s | p.T22T22z;
    s | p.T22T22z_sus;

    s | p.T20T20x;
    s | p.T20T20x_sus;
    s | p.T22T22x;
    s | p.T22T22x_sus;

    s | p.T20T20y;
    s | p.T20T20y_sus;
    s | p.T22T22y;
    s | p.T22T22y_sus;

    s | p.T32T32;
    s | p.T32T32_sus;

    s | p.mean_autocorrelation_time;

    s | p.kappa;
}

#endif  /* _STANDARDPROPERTIES_H */


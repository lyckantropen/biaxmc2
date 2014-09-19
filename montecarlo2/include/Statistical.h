/*
 * File:   Statistical.h
 * Author: karol
 *
 * Created on 19 listopad 2009, 15:02
 *
 * Several types and algorithms for calculating
 * averages, errors and fluctuations.
 *
 *
 */



#ifndef _STATISTICAL_H
#define _STATISTICAL_H

#include <iostream>
#include <cmath>
#include "Maths.h"
#include "serializer.h"
#include "Random.h"


/**
 * @brief This type defines a real variable with a definite error.
 * All arithmetic is transparent and identical to that of a 'double'
 * variable and acts directly upon the (expectation) value. The error
 * can be accessed via the Error() function. The value of the standard
 * error is calculated at each calculation automatically with respect
 * to the gradient principle.
 *
 */
class Value
{
    ///define the serialization operator
    template<class serializer_t>
    friend void operator|(serializer_t &, Value&);

    double val;
    double err;
public:
    Value();
    Value(const double & v, const double & e);
    Value(const double & v);

    operator double & ();
    operator const double & () const;
    const Value operator=(const double & v);
    const Value operator+(const Value & v) const;
    const Value operator-(const Value & v) const;
    const Value operator*(const Value & v) const;
    const Value operator/(const Value & v) const;
    bool operator==(const double & v);
    bool operator==(const Value & v);

    ///get the error
    const double & Error() const;
    ///print the value and error as a formatted string
    const std::string Print() const;
    ///print the value and error separated by a tab
    const std::string TableForm() const;
    ///print the value and error as Mathematica array
    const std::string MathematicaForm() const;
};

///serialization operator
template<class serializer_t>
void operator|(serializer_t & s, Value & v)
{
    s | v.val;
    s | v.err;
}

///this square root function automatically calculates the error
inline Value sqrt(const Value & v)
{
    double sq = std::sqrt(std::abs(double(v)));
    return Value(sq, v.Error() / 2 / sq);
}

///the output stream operator prints v as a formatted string
extern std::ostream & operator<<(std::ostream & s, Value & v);

/*
template<class value_type>
Value MeanSub(const std::valarray<value_type> & v){
    int limit = v.size();
    int start = 0;
    value_type mean = v.sum()/value_type(v.size());
    std::valarray<value_type> dev = std::pow(v-mean,2);
    return Value(mean,std::sqrt(dev.sum()/(value_type(v.size())-value_type(1.0)) ));
}
 */

/// Calculates the arithmetic average and returns it as a Value with the error
/// calculated as a fluctuation over a fixed number of 'slices': the sample 'v' is
/// divided into '_slices' ranges (if larger) and '_slices' averages are calculated
/// from them. The fluctuation is calculated as a square root of variance of these
/// partial averages.
template<class value_type>
Value Mean(const std::valarray<value_type> & v, const int start = 0, const int limit = 0, const int _slices = 100)
{
    int lim;
    int slices = _slices;
    if(limit == 0)
        lim = v.size();
    else
        lim = limit;
    if(slices > v.size())
    {
        slices = v.size();
        //    std::cout << v.size() << ", but slices: " << _slices << ". Decimated to " << slices << std::endl;
    }
    if(slices == 0) slices = 1;
    /*
    std::valarray<value_type> relevant = v[std::slice(start,lim-start,1)];
    value_type mean = relevant.sum()/value_type(relevant.size());
    std::valarray<value_type> dev = std::pow(relevant-mean,2);
    return Value(mean,std::sqrt(dev.sum()/(value_type(relevant.size())-value_type(1.0)) ));
     */
    std::valarray<value_type> means(0.0, slices);
    int samplesize = (lim - start) / slices;
    for(int i = 0; i < slices; i++)
    {
        std::valarray<value_type> submean = v[std::slice(start + i * samplesize, samplesize, 1)];
        means[i] = submean.sum() / value_type(samplesize);
    }
    value_type mean = means.sum() / value_type(slices);
    std::valarray<value_type> dev = std::pow(means - mean, 2);
    return Value(mean, std::sqrt(dev.sum() / value_type(slices)));
}

/// Calculates an average array from an array of arrays of identical dimension.
/// Note: does not calculate standard deviation.
template<class value_type>
std::valarray<value_type> MeanVector(const std::vector<std::valarray<value_type> > & v, const int start = 0, const int limit = 0)
{
    int minn, maxn;
    if(limit == 0)
        maxn = v.size();
    else
        maxn = limit;
    std::valarray<value_type> result(0.0, v[0].size());
    for(int i = start; i < maxn; i++)
    {
        //    std::cout << "v["<<i<<"]=" << v[i] << std::endl;
        result += v[i];
    }
    result /= double(maxn - start);
    //std::cout << "MeanVector: " << result << std::endl;

    return result;
}

/// A bootstrap variant of the Mean function above. The fluctuation is calculated
/// over 'res' random samples of the provided sample 'v'. Although numerous
/// checks have been made, I find that the calculated fluctuation is dubiously small
/// and decreases with sample size. As of now it is not used throughout biaxmc2,
/// so please verify before use.
//Efron, Tibshirani, "An Introduction to Bootstrap", 1993 Chapman & Hall, Inc., pp. 45-47
template<class value_type>
Value BootstrapMean(const std::valarray<value_type> & v, const int start = 0, const int limit = 0, const int res = 200)
{
    int resamples = res;
    if(res == 0)
        resamples = v.size();
    int lim;
    if(limit == 0)
        lim = v.size();
    else lim = limit;
    // wybieramy zakres zmiennej, który nas interesuje
    std::valarray<value_type>   relevant = v[std::slice(start, lim - start, 1)];
    int N = relevant.size();
    // estymatory wartości oczekiwanej liczone dla każdej próbki
    std::valarray<value_type> bootmean(resamples);
    value_type B = value_type(resamples);

    for(int i = 0; i < resamples; i++)
    {
        value_type * sample = new value_type[N];
        // próbkowanie z powtórzeniami
        for(int j = 0; j < N; j++)
        {
            sample[j] = relevant[int(N * random01())];
        }
        bootmean[i] = std::valarray<value_type>(sample, N).sum() / value_type(N);
        delete sample;
    }
    value_type Bsum = bootmean.sum() / B;
    double stddev = std::sqrt(std::pow((bootmean - Bsum), value_type(2.0)).sum() / (B - value_type(1)));
    return Value(bootmean.sum() / B, stddev);
}

/// Calculates fluctuation of the array 'variable' up to the index 'acc_idx'
/// with respect to the multiple sampled average algorithm, where the number
/// of resampling iterations is equal to the size of the original sample.
extern Value CalculateFluctuation(const vect & variable, const int & acc_idx = 0);

//#if defined(__GNUC__) && (__GNUC_MINOR__ < 7)
//template Value Mean<Value>(const std::valarray<Value> & v, const int, const int, const int);
//template  Value Mean<double>(const std::valarray<double> & v, const int, const int, const int);
//template  std::valarray<double> MeanVector<double>(const std::vector<vect> & v, const int start, const int limit);
//template  Value BootstrapMean<double>(const std::valarray<double> & v, const int start, const int limit, const int resamples);
//template<>  std::valarray<double> MeanVector<>(const std::vector<std::valarray<Value> > & v, const int & limit);
//#endif //__APPLE__
#endif  /* _STATISTICAL_H */


/* 
 * File:   Statistical.h
 * Author: karol
 *
 * Created on 19 listopad 2009, 15:02
 */

#ifndef _STATISTICAL_H
#define	_STATISTICAL_H

#include <iostream>
#include <cmath>
#include "valarray_external.h"
#include "serializer.h"
#include "random01.h"

class Value {
    template<class serializer_t>
    friend void operator|(serializer_t &,Value&);
    double val;
    double err;
public:
    Value():val(0),err(0) {}
    Value(const double & v, const double & e){
        val=v;
        err=e;
    }
    Value(const double & v){
        val=v;
        err=0.0;
    }
    
    operator double & (){
        return val;
    }
    operator const double & () const {
        return val;
    } 
    const double & Error() const {
        return err;
    }
    const Value operator=(const double & v){
        val=v;
        return *this;
    }
    const Value operator+(const Value & v) const {
        return Value(val+v.val,err+v.err);
    }
    const Value operator-(const Value & v) const {
        return Value(val-v.val,err+v.err);
    }
    const Value operator*(const Value & v) const {
        return Value(val*v.val,std::abs(err*v.val)+std::abs(v.err*val));
    }
    const Value operator/(const Value & v) const {
        return Value(val/v.val,std::abs(err/v.val)+std::abs(v.err/val/val));
    }
    bool operator==(const double & v){
        return val==v;
    }
    bool operator==(const Value & v){
        return val==v.val;
    }
    const std::string Print() const {
        std::stringstream s;
        s << std::setprecision(2);
        s << val << "("<<err<<","<< std::abs(err/val*100) << "%)";
        return s.str();
    }
    const std::string TableForm() const {
        std::stringstream s;
        s << val << "\t" << err ;
        return s.str();
    }
};

//extern void operator|(boostbase::serializer &,Value&);
template<class serializer_t>
void operator|(serializer_t & s,Value & v){
    s|v.val;
    s|v.err;
}
inline Value sqrt(const Value & v){
    double sq = std::sqrt(std::abs(double(v)));
    return Value(sq,v.Error()/2/sq);
}
extern std::ostream & operator<<(std::ostream & s, Value & v);


template<class value_type>
Value Mean(const std::valarray<value_type> & v,const int start=0, const int limit=0){
    int lim;
    if(limit==0)
        lim=v.size();
    else
        lim=limit;
    std::valarray<value_type> relevant = v[std::slice(start,lim-start,1)];
    value_type mean = relevant.sum()/value_type(relevant.size());
    std::valarray<value_type> dev = std::pow(relevant-mean,2);
    return Value(mean,std::sqrt(dev.sum()/(value_type(relevant.size())-value_type(1.0)) ));
}

template<class value_type>
std::valarray<value_type> MeanVector(const std::vector<std::valarray<value_type> > & v,const int start=0, const int limit=0){
    int minn,maxn;
    if(limit==0)
        maxn=v.size();
    else
        maxn=limit;
    std::valarray<value_type> result(v[0].size());
    for(int i=start;i<maxn;i++){
//        std::cout << "v["<<i<<"]=" << v[i] << std::endl;
        result+=v[i];
    }
    result/=double(maxn);
//    std::cout << "MeanVector: " << result << std::endl;

    return result;
}
//Efron, Tibshirani, "An Introduction to Bootstrap", 1993 Chapman & Hall, Inc., pp. 45-47
template<class value_type>
Value BootstrapMean(const std::valarray<value_type> & v, const int start=0, const int limit=0, const int resamples=50){
    int lim;
    if(limit==0)
        lim=v.size();
    else lim=limit;
    // wybieramy zakres zmiennej, który nas interesuje
    std::valarray<value_type>   relevant = v[std::slice(start,lim-start,1)];
    int N = relevant.size();
    // estymatory wartości oczekiwanej liczone dla każdej próbki
    std::valarray<value_type> bootmean(resamples);
    value_type B = value_type(resamples);

    for(int i=0;i<resamples;i++){
        value_type * sample = new value_type[N];
        // próbkowanie z powtórzeniami
        for(int j=0;j<N;j++){
            sample[j]=relevant[int(N*random01())];
        }
        bootmean[i]=std::valarray<value_type>(sample,N).sum()/value_type(N);
        delete sample;
    }
    value_type Bsum = bootmean.sum()/B;
    double stddev = std::sqrt(std::pow((bootmean-Bsum),value_type(2.0)).sum()/(B-value_type(1)));
    return Value(bootmean.sum()/B,stddev);
};


template  Value Mean<Value>(const std::valarray<Value> & v,const int, const int);
template  Value Mean<double>(const std::valarray<double> & v, const int, const int);
template  std::valarray<double> MeanVector<double>(const std::vector<vect> & v, const int start, const int limit);
template  Value BootstrapMean<double>(const std::valarray<double> & v, const int start, const int limit, const int resamples);
//template<>  std::valarray<double> MeanVector<>(const std::vector<std::valarray<Value> > & v, const int & limit);

#endif	/* _STATISTICAL_H */


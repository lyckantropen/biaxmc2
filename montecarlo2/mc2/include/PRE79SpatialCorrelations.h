/* 
 * File:   PRE79SpatialCorrelations.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 12:00
 */

#ifndef _PRE79SPATIALCORRELATIONS_H
#define	_PRE79SPATIALCORRELATIONS_H

#include "SpatialCorrelationEvolution.h"

class Delta200Correlation:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta200Correlation &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 3.0/2.0*Rank2Contraction((a.GetQZ()-Identity(3)/3.0),(b.GetQZ()-Identity(3)/3.0));
    }
public:
    Delta200Correlation(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};


class Delta220Correlation:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta220Correlation &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return std::sqrt(3.0)/2.0*Rank2Contraction((a.GetQX()-a.GetQY()),(b.GetQZ()-Identity(3)/3.0));
    }
public:
    Delta220Correlation(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta222Correlation:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta222Correlation &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 0.5*Rank2Contraction((a.GetQX()-a.GetQY()),(b.GetQX()-b.GetQY()));
    }
public:
    Delta222Correlation(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};


class Delta322Correlation:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta322Correlation &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return Rank3Contraction(a.GetT(),b.GetT());
    }
public:
    Delta322Correlation(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};


class ParityCorrelation:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, ParityCorrelation &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return a.GetParity()*b.GetParity();
    }
public:
    ParityCorrelation(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

template<class serializer_t>
void operator|(serializer_t & s, Delta200Correlation & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta220Correlation & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta222Correlation & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta322Correlation & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, ParityCorrelation & e){
    e.SerializerHelper<serializer_t>(s);
}



#endif	/* _PRE79SPATIALCORRELATIONS_H */


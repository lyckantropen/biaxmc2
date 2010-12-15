/* 
 * File:   PRE79SpatialCorrelations.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 12:00
 */

#ifndef _PRE79SPATIALCORRELATIONS_H
#define	_PRE79SPATIALCORRELATIONS_H

#include "SpatialCorrelationEvolution.h"

/*
 * Funkcje G 2,0,0 względem osi Z,X i Y
 */

class Delta200CorrelationZ:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta200CorrelationZ &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 3.0/2.0*Rank2Contraction((a.GetQZ()-Identity(3)/3.0),(b.GetQZ()-Identity(3)/3.0));
    }
public:
    Delta200CorrelationZ(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta200CorrelationX:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta200CorrelationX &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 3.0/2.0*Rank2Contraction((a.GetQX()-Identity(3)/3.0),(b.GetQX()-Identity(3)/3.0));
    }
public:
    Delta200CorrelationX(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta200CorrelationY:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta200CorrelationY &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 3.0/2.0*Rank2Contraction((a.GetQY()-Identity(3)/3.0),(b.GetQY()-Identity(3)/3.0));
    }
public:
    Delta200CorrelationY(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

/*
 * Funkcje G 2,2,0 względem osi Z,X,Y
 */

class Delta220CorrelationZ:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta220CorrelationZ &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return std::sqrt(3.0)/2.0*Rank2Contraction((a.GetQX()-a.GetQY()),(b.GetQZ()-Identity(3)/3.0));
    }
public:
    Delta220CorrelationZ(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta220CorrelationX:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta220CorrelationX &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return std::sqrt(3.0)/2.0*Rank2Contraction((a.GetQY()-a.GetQZ()),(b.GetQX()-Identity(3)/3.0));
    }
public:
    Delta220CorrelationX(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta220CorrelationY:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta220CorrelationY &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return std::sqrt(3.0)/2.0*Rank2Contraction((a.GetQZ()-a.GetQX()),(b.GetQY()-Identity(3)/3.0));
    }
public:
    Delta220CorrelationY(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

/*
 * Funkcje G 2,2,2 względem osi Z,X,Y
 */

class Delta222CorrelationZ:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta222CorrelationZ &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 0.5*Rank2Contraction((a.GetQX()-a.GetQY()),(b.GetQX()-b.GetQY()));
    }
public:
    Delta222CorrelationZ(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta222CorrelationX:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta222CorrelationX &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 0.5*Rank2Contraction((a.GetQY()-a.GetQZ()),(b.GetQY()-b.GetQZ()));
    }
public:
    Delta222CorrelationX(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

class Delta222CorrelationY:public SpatialCorrelationEvolution {
    using SpatialCorrelationEvolution::SerializerHelper;
    template<class serializer_t>
    friend void operator|(serializer_t &, Delta222CorrelationY &);
    virtual double CalculateContraction(const Particle &a, const Particle &b){
        return 0.5*Rank2Contraction((a.GetQZ()-a.GetQX()),(b.GetQZ()-b.GetQX()));
    }
public:
    Delta222CorrelationY(Lattice * l=NULL, int nc=0):SpatialCorrelationEvolution(l,nc){}
};

/**
 * Funkcja korelacji G 3,2,2
 */
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

/*
 * Serializery dla funkcji G 2,0,0
 */
template<class serializer_t>
void operator|(serializer_t & s, Delta200CorrelationZ & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta200CorrelationX & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta200CorrelationY & e){
    e.SerializerHelper<serializer_t>(s);
}
/*
 * Serializery dla funkcji G 2,2,0
 */
template<class serializer_t>
void operator|(serializer_t & s, Delta220CorrelationZ & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta220CorrelationX & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta220CorrelationY & e){
    e.SerializerHelper<serializer_t>(s);
}

/*
 * Serializer dla funkcji G 2,2,2
 */
template<class serializer_t>
void operator|(serializer_t & s, Delta222CorrelationZ & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta222CorrelationX & e){
    e.SerializerHelper<serializer_t>(s);
}
template<class serializer_t>
void operator|(serializer_t & s, Delta222CorrelationY & e){
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


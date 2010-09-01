#include "valarray_external.h"
#include <cmath>

std::ostream & operator<<(std::ostream & o, const vect & v){
    for(int i=0;i<v.size();i++){
        o << v[i] << " ";
    }
    return o;
}

std::string MathematicaForm(const vect & v){
    std::stringstream s;
    s << std::setprecision(15) << std::fixed;
    s << "{";
    for(int i=0;i<v.size();i++){
        s << v[i] ;
        if(i<v.size()-1)
            s << ",";
    }
    s<< "}";
    return s.str();
}

double  DotProduct(const vect & a, const vect & b){
    return (a*b).sum();
}

double MatrixDotProduct(const vect & a, const vect & b){
    vect p(0.0,6);
    p = a*b;
    return p[0]+p[3]+p[5]+2.*p[1]+2.*p[2]+2.*p[4];
}

double  Norm(const vect & v){
    vect result = std::pow(v,2.0);
    return std::sqrt(result.sum());
}
vect Identity(const int & dim=3){
    //TODO: wyższe wymiary niż 3
    vect out(6);
    out[0]=out[3]=out[5]=1.0;
    out[1]=out[2]=out[4]=0.0;
    return out;
}

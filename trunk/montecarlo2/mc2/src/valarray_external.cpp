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
    return a[0]*b[0]+
	   a[3]*b[3]+
	   a[5]*b[5]+
	   2.0*a[1]*b[1]+
	   2.0*a[2]*b[2]+
	   2.0*a[4]*b[4];
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

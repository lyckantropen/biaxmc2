#include "Maths.h"
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








double Minimum(const vect & a){
    double min=a[0];
    for(int i=0;i<a.size();i++){
        if(a[i]<min) min=a[i];
    }
    return min;
}

int MinimumIndex(const vect & a){
    double min=a[0];
    int mid=0;
    for(int i=0;i<a.size();i++){
        if(a[i]<min){
            min=a[i];
            mid=i;
        }
    }
    return mid;
}

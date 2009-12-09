#include "Particle.h"

std::ostream & operator<<(std::ostream & s, const Particle & p){
    s<<p.parity ;
    for(int i=0;i<p.x.size();i++){
        s << " "<< p.x[i] ;
    }
    return s;
}



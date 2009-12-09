/* 
 * File:   Lattice.cpp
 * Author: karol
 * 
 * Created on 17 listopad 2009, 17:41
 */

#include "Lattice.h"
std::ostream & operator<<(std::ostream & s,const Lattice & lat){
    for(int i=0;i<lat.Particles.size();i++){
        s << lat.Particles[i] << std::endl;
    }
    s << std::endl;
    return s;
}


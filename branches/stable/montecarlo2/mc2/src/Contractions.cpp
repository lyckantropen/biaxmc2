#include "Contractions.h"

double  Rank2Contraction(const vect & a, const vect & b){
        vect coeff(6);
        coeff[0]=coeff[3]=coeff[5]=1.0;
        coeff[1]=coeff[2]=coeff[4]=2.0;
        coeff*=a*b;
        return coeff.sum();
}

// kantrakcja tensora 3x3x kodowanego 10-cioma składowymi
double  Rank3Contraction(const vect & a, const vect & b){
        vect coeff(10);
        coeff[0]=1;
        coeff[1]=3;
        coeff[2]=3;
        coeff[3]=1;
        coeff[4]=3;
        coeff[5]=6;
        coeff[6]=3;
        coeff[7]=3;
        coeff[8]=3;
        coeff[9]=1;
        coeff*=a*b;
        return coeff.sum();
}


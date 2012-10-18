#include "Statistical.h"

Value CalculateFluctuation(const vect & variable,const int & acc_idx){
    int size=acc_idx+1;
    if(acc_idx==0)
        size=variable.size();

    vect fluct3(size);
    for(int i=0;i<(size);i++){
        vect e(size);
        vect e2(size);
        for(int j=0;j<(size);j++){
            int t = (size)*random01();
            e[j]=variable[t];
            e2[j]=variable[t]*variable[t];
        }
        fluct3[i]=(double(Mean(e2))-double(Mean(e))*double(Mean(e)));
    }
    return Mean(fluct3);
}

std::ostream & operator<<(std::ostream & s, Value & v){
    s << double(v) << " " << v.Error() ;
    return s;
}


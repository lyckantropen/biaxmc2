#include "Statistical.h"

Value CalculateFluctuation(const vect & variable,const int & acc_idx){
        vect fluct3(acc_idx+1);
        for(int i=0;i<(acc_idx+1);i++){
            vect e(acc_idx+1);
            vect e2(acc_idx+1);
            for(int j=0;j<(acc_idx+1);j++){
                int t = (acc_idx)*random01();
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


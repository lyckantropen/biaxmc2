#include "4DSphereRW.h"



vect    RandomWalkOn4DSphere(const double & r,const vect & x){
    vect step = RandomPointOn4DSphereMarsaglia(r);
    //vect step = RandomPointOn4DSphereOld(r);
    vect result = x + step;
    return result/Norm(result);
}

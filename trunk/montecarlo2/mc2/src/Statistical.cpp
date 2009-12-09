#include "Statistical.h"

std::ostream & operator<<(std::ostream & s, Value & v){
    s << double(v) << " " << v.Error() ;
    return s;
}


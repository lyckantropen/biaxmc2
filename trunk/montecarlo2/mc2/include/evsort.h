#ifndef _EVSORT_H_

#include "std.h"
#include "boost.h"
#include "valarray_external.h"

struct evs {
	double e;
	double v[3];
	evs() {}
	evs(const double _e, const double _v[3]){
		e=_e;
		v[0]=_v[0];
		v[1]=_v[1];
		v[2]=_v[2];
	}
};

inline bool operator<(const evs & a, const evs &b){
	return std::abs(a.e)>std::abs(b.e);
}


#endif

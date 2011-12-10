/* 
 * File:   random01.cpp
 * Author: karol
 * 
 * Created on 17 listopad 2009, 12:50
 */

#include "random01.h"


boost::rand48   rng2;
//rangen_wd random01;
//rangen random01;
//MarsagliaRNG random01;
//rangen_mwc random01;
randtype random01;
template<class rng> std::vector<rng> rng_wrap<rng>::rg;
//rangen * rg = NULL;

/*
double random01(){
    double retval = 0.0;
    if(rg==NULL)
        rg = new rangen;

    retval = rg->Gen01();
    return retval;
}
*/
int plusminusone(){
    if(random01()>0.5)
        return 1;
    else return -1;
}
///this algorithm is due to Marsaglia
vect    RandomPointOn4DSphereMarsaglia(const double & r){
    vect result(4);
    double r1=0,r2=0;
    double y1=0,y2=0,y3=0,y4=0;
    do {
        y1=(1.0-2.0*random01());
        y2=(1.0-2.0*random01());
        r1=y1*y1+y2*y2;
    } while(r1>1);
    do {
        y3=(1.0-2.0*random01());
        y4=(1.0-2.0*random01());
        r2=y3*y3+y4*y4;
    } while(r2>1);
    /*
     * This is OK, in this notation r1 and r2 are squares of length, as required
     */
    double sr=std::sqrt((1-r1)/r2);

    result[0]=r*y1;
    result[1]=r*y2;
    result[2]=r*y3*sr;
    result[3]=r*y4*sr;
    return result;
}

vect    RandomPointOn4DSphereOld(const double & r){
	vect x(4);
	double r1;
	do
	{
		x[0]=(1-2*random01());
		x[1]=(1-2*random01());
		x[2]=(1-2*random01());
		x[3]=(1-2*random01());
		r1=Norm(x);
	} while (r1>1.0);
	x[0]=x[0]/r1*r;
	x[1]=x[1]/r1*r;
	x[2]=x[2]/r1*r;
	x[3]=x[3]/r1*r;
	return x;
}

long int seed1 = 46723;			// seed for random long integer generator
long int seed2 = 2593459;		// seed for shuffling algorithm

long int jj =89343;
double denom=(1.0/(m-1));
const int NR=64;

long int y;			// random int mixed into shuffling array.
long int j[NR];		// vector of random integers that are shuffled
					// Both y and j[] are seeded by making the singl

long int random_int(long int& j)
{
	long int l;
	l = j/q;				// Parameter "j" is recursively defined to produce a
							// a sequence of random numbers.
	j = a*(j - q*l) -r*l;	// This iteration step is actually a mod 2^31
							// operation as described by Schrage.  Here, first term =
							// j mod q and second term is (int divide(j/q) )*r.
							// See Barkema, page 388.
	if (j<0) j += m;		// This shift bring the number back into the interval
							// [1, m - 1].
	return (j-1);	// j - 1 is a random integer in interval [0, m - 2]
}	//end program to return random long integer

double drandom(long int& seed2)
{
	long int l;
	long int k;

	l = seed2/q;
	seed2 = a*(seed2 - q*l) - r*l;	//	Selects "seed2" to be a random integer in
							// [1, 2^31 - 2 ] and updates "seed2" as reference variable.
							//	This ensures that a series of calls will produce a
							//	a series of random long ints "seed2".
	if (seed2 < 0)			//ensures that seed2 is random int in [1, m - 1]
		{
			seed2 = seed2 + m;
		}
	k = (y/m)*NR;		//Selects an integer in [0, NR]

	y = j[k];			// Sets y = integer of k-th component of j[].
						// Side effect: global variable y is updated for
						// next call of drandom.
	j[k] = seed2;		// Side effect: global variable j[k] is updated
						// for next call of drandom.

	return denom*(y - 1);	// Returns randomly shuffled double float in [0, 1).
}	//end Shuffler program

void initializeRNG(long int& seed1)
{
	int s;
	for(s = 0; s< NR ; s++)	//loop iteratively seeds each component of j[].
		{
			j[s] = random_int(seed1); //seeds component "s" of shuffler vector j[].
			//	std::cout << "j[" << s << "] = " << j[s] << '\n';
			//	The above line is used as a test to insure that j[] is a random array
			//	That is, to print out array j[] to check random initialization
		}
	y = random_int(seed1);	// seeds variable j
	//	std::cout << "yy = " << yy<< '\n';
	//	That is, to print out yy to check random initialization
}	//end program to initialize Shuffler program


//template<> rangen * Singleton<rangen>::instance=NULL;

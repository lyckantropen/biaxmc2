/* 
 * File:   WydroRNG.h
 * Author: karol
 *
 * Created on 27 styczeń 2011, 11:11
 */

#ifndef _WYDRORNG_H
#define	_WYDRORNG_H

//Declarations and initializations for random integer generator

const long int m=2147483647; // m-1 = number of random integers generable
const long int a=16807;		 // values of constants a, q, and r give good
const long int q=127773;	 // randomization and maximum cyle length m-1.
const long int r=2836;
//Additional declarations and initializations needed for shuffler program
extern long int seed1;			// seed for random long integer generator
extern long int seed2;		// seed for shuffling algorithm

//minstd_rand0 16807, 0, 2147483647, 1043618065
//minstd_rand 48271, 0, 2147483647, 399268537



//_____________________________________________________________________________________________

extern long int random_int(long int& j);



/* B) Remainder for Shuffler program that generates double float in [0,1). */

// PRE-COND:	Variable seed2 has been declared and initialized in [1, 2^31 - 2].
//				Also, a single call must be made to initializeRNG(seed1) to seed the
//				y and j[] variables of the Shuffler program.

// POST-COND:	Shuffler returns random double float in [0, 1).  Repeat cycle is
//				believed to be about (m - 1)^64 where m = 2^31 - 1.  Thus, m is
//				about = 2x10^9, and the repeat cycle, i.e., (m - 1)^64 is thus,
//				more or less infinite in any reasonably doable simulation.
//				Including this file in a program makes each call of
//				drandom(long int& i) update i, yy, and j[k] as a side effect.  The
//				reason is: seed2, y, and j[] are global variables.

//Shuffler program uses the generator of random integers and declarations therein.


/*
extern long int jj;
extern double denom;
const int NR = 64;

extern long int y;			// random int mixed into shuffling array.
extern long int j[];		// vector of random integers that are shuffled
					// Both y and j[] are seeded by making the single call
					// of initializeRNG(seed1).

*/


//_____________________________________________________________________________________________

extern double drandom(long int& seed2);

/* C) Program to initialize Shuffler program. */

// PRE-COND:	Variable seed1 has been declared and initialized in [1, 2^31 - 2].
//				Note: variable y and vector j[] are globally declared in a program
//				that includes this header file.
// POST-COND:	Shuffler variable y and components of vector j[] are seeded with
//				randomly selected long integers in [1, m - 1].
//				Side effect of fact that seed1 is reference variable is that
//				each initialization of a component of j[] updates seed1 thereby
//				ensuring that the sequence of components are given a random sequence
//				of long int values.






//_____________________________________________________________________________________________

extern void initializeRNG(long int& seed1);

#endif	/* _WYDRORNG_H */


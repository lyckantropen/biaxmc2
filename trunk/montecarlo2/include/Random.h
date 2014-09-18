/*
 * File:   Random.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 12:50
 */

#ifndef _RANDOM01_H
#define _RANDOM01_H
#include "boost.h"
#include "std.h"
#include "WydroRNG.h"
#include "Maths.h"
#include "randomc.h"    ///multiply with carry

/// Mersenne Twister 19937 wrapper.
/// The RNG is seeded from a rand48 RNG for each thread,
/// however this does not produce satisfyingly uncorrelated
/// sequences.
///
/// Not used.
///
class rangen
{
    boost::mt19937  rng;
    boost::uniform_01<boost::mt19937>   uni01;

public:
    rangen();
    /// The copy-constructor is launched at the OpenMP fork
    /// for each thread when we declare that random01 should
    /// be private. However, we need the individual RNGs
    /// to be initialized again, not copied.
    ///
    rangen(const rangen &);
    double Gen01();             ///compat.
    double operator()();        ///returns random real on [0,1)
};
extern boost::rand48 rng2;

/// This one uses the linear congruential generator from
/// Tomasz Wydro
///
/// Not thread-independent. Not used.
///
class rangen_wd
{
public:
    rangen_wd();
    double Gen01();
    double operator()();
};

/// A thread-safe OpenMP-based wrapper for the Multiply With Carry
/// RNG by George Marsaglia.
/// The thread limit is 128.
///
/// Just use a separate RNG for each OpenMP thread.
///
class rangen_mwc
{
    CRandomMother gen;
    unsigned long lsp[128][3];
public:
    rangen_mwc();
    /// The copy-constructor is launched at the OpenMP fork
    /// for each thread when we declare that random01 should
    /// be private. However, we need the individual RNGs
    /// to be initialized again, not copied.
    ///
    rangen_mwc(const rangen_mwc &);
    void setup(int k);
    double operator()();     ///returns random real on [0,1)
};

/// This is a wrapper for the 'rangen_mwc' RNG
/// which explicitly allocates 'n' RNGs for 'n' threads.
/// Apparently is slow. Not used.
class mwc_wrap
{
public:
    static std::vector<rangen_mwc> rg;
    static void setup(int n);
    double operator()();

};

/// If you want to use a different RNG, uncomment the desired one.
//extern rangen_wd random01;
//extern rangen random01;
extern rangen_mwc random01;
//extern mwc_wrap random01;

/// This just returns 1 or -1 with the same probability
int plusminusone();


#endif  /* _RANDOM01_H */


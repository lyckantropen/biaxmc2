#include "Metropolis.h"
#include "Hamiltonian.h"
#include "4DSphereRW.h"
#include "Lattice.h"

Metropolis::Metropolis(const Settings &set, const Hamiltonian * h):
    settings(&set),
    hamiltonian(h),
    radius(set.simulation.radius),
    acc_llimit(set.simulation.metropolis_lower_acceptance_limit),
    acc_ulimit(set.simulation.metropolis_higher_acceptance_limit),
    parity_prob(set.simulation.parity_flip_probability)
{

}

Metropolis::Metropolis(const Metropolis &s):
    settings(s.settings)
{
    hamiltonian = s.hamiltonian;
    radius = s.radius;
    acc_llimit = s.acc_llimit;
    acc_ulimit = s.acc_ulimit;
    parity_prob = s.parity_prob;
}

const Metropolis &Metropolis::operator=(const Metropolis &s)
{
    hamiltonian = s.hamiltonian;
    radius = s.radius;
    acc_llimit = s.acc_llimit;
    acc_ulimit = s.acc_ulimit;
    parity_prob = s.parity_prob;
    settings = s.settings;
    return *this;
}

vect Metropolis::OrientationNudge(const vect &old) const
{
    return RandomWalkOn4DSphere(radius, old);
}

short Metropolis::ParityNudge(const short &old) const
{

    if(random01() < parity_prob)
        return -old;
    else return old;

    //return plusminusone();
}

bool Metropolis::Accept(const double &dE) const
{
    if(dE < 0)
        return true;
    else
    {
        double acceptance = std::exp(-dE);
        if(random01() < acceptance)
            return true;
        else
            return false;
    }

}

const Hamiltonian *Metropolis::GetHamiltonian() const
{
    return hamiltonian;
}

double Metropolis::MeasureAccepted(Lattice *lat)
{
    if(lat == NULL) return -1 ;
    double N = 0.0;
    int acc_moves = 0;
    Lattice testlat = *lat;

    double tries = 5;
    N = testlat.GetN();
    for(int i = 0; i < tries; i++)
    {
        int acc_rot = 0, acc_p = 0;
        testlat.Sweep(this, acc_rot, acc_p);
        acc_moves += acc_p;
    }
    return double(acc_moves) / (N * tries);

}

void Metropolis::AdjustRadius(Lattice *lat, const double &decimation)
{
    if(!settings->simulation.adjust_radius) return;

    if(lat == NULL) return ;
    double acc_fraction = 0.0;
    double N = 0.0;
    int acc_moves = 0;


    while(acc_fraction < acc_llimit || acc_fraction > acc_ulimit)
    {
        Lattice testlat = *lat;
        acc_moves = 0;
        double tries = 2;
        N = testlat.GetN();
        for(int i = 0; i < tries; i++)
        {
            int acc_rot = 0, acc_p = 0;
            testlat.Sweep(this, acc_rot, acc_p);
            acc_moves += acc_rot;
        }
        acc_fraction = double(acc_moves) / (N * tries);


        if(acc_fraction < acc_llimit)
        {
            radius *= (1.0 - decimation);
            //Log() << "decimating down to r=" << radius << " because acc_fraction=" << acc_fraction << std::endl;
        }
        if(acc_fraction > acc_ulimit)
        {
            radius *= (1.0 + decimation);
            //Log() << "decimating up to r=" << radius << " because acc_fraction=" << acc_fraction << std::endl;
        }

        if(radius >= 1.0)
        {
            radius = 0.999;
            break;
        }
        if(radius <= 0.001)
        {
            radius = 0.001;
            break;
        }

    }
    Log() << acc_fraction << " (" << acc_moves << "/" << N << ") accepted, radius " << radius << std::endl;
}

void Metropolis::SetStream(std::ostream *os)
{
    ILoggable::SetStream(os);
}

const double &Metropolis::GetAccLLimit() const
{
    return acc_llimit;
}

const double &Metropolis::GetAccULimit() const
{
    return acc_ulimit;
}

const double &Metropolis::GetParityProb() const
{
    return parity_prob;
}

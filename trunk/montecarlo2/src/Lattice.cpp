/*
 * File:   Lattice.cpp
 * Author: karol
 *
 * Created on 17 listopad 2009, 17:41
 */

#include "Lattice.h"
std::ostream & operator<<(std::ostream & s, const Lattice & lat)
{
    for(int i = 0; i < lat.Particles.size(); i++)
    {
        s << lat.Particles[i] << std::endl;
    }
    s << std::endl;
    return s;
}
void Lattice::Construct()
{
    Log() << "Constructing\n";
    if(Particles.size() != 0) Particles.clear();
    Particles.resize(N);

    for(int i = 0; i < N; i++)
    {
        Particle & cur = Particles[i];
        int A = L * W;
        int x = 0, y = 0, z = 0;
        x = i % L;
        y = i % A / L;
        z = i / A;
        cur.PlaceAt(x, y, z);

        if(x == L - 1 && periodic_L)
        {
            cur.Connect(Particles[i - L + 1], i - L + 1);
            Log() << "L boundary: " << i << " at " << cur.GetR() << " to " << i - L + 1 << std::endl;
        }
        else if(i + 1 < N) cur.Connect(Particles[i + 1], i + 1);

        if(x == 0 && periodic_L)
        {
            cur.Connect(Particles[i + L - 1], i + L - 1);
            Log() << "L boundary: " << i << " at " << cur.GetR() << " to " << i + L - 1 << std::endl;
        }
        else if(i - 1 >= 0) cur.Connect(Particles[i - 1], i - 1);

        if(y == W - 1 && periodic_W)
        {
            cur.Connect(Particles[i - A + L], i - A + L);
            Log() << "W boundary: " << i << " at " << cur.GetR() << " to " << i - A + L << std::endl;
        }
        else if(i + L < N) cur.Connect(Particles[i + L], i + L);

        if(y == 0 && periodic_W)
        {
            cur.Connect(Particles[i + A - L], i + A - L);
            Log() << "W boundary: " << i << " at " << cur.GetR() << " to " << i + A - L << std::endl;
        }
        else if(i - L >= 0) cur.Connect(Particles[i - L], i - L);

        if(z == H - 1 && periodic_H)
        {
            cur.Connect(Particles[i - N + A], i - N + A);
            Log() << "H boundary: " << i << " at " << cur.GetR() << " to " << i - N + A << std::endl;
        }
        else if(i + A < N) cur.Connect(Particles[i + A], i + A);

        if(z == 0 && periodic_H)
        {
            cur.Connect(Particles[i + N - A], i + N - A);
            Log() << "H boundary: " << i << " at " << cur.GetR() << " to " << i + N - A << std::endl;
        }
        else if(i - A >= 0) cur.Connect(Particles[i - A], i - A);

    }
    for(int i = 0; i < N; i++)
    {
        Log() << "Particle " << i << " at " << Particles[i].GetR() << " has neighbors: \n";
        for(int n = 0; n < Particles[i].GetNNeighbors(); n++)
        {
            Log() << Particles[i].GetNeighborIndex(n) << ": " << Particles[i].GetNeighbor(n)->GetR() << ", " ;
        }
        Log() << std::endl;
    }

}

Lattice::Lattice(): L(0), W(0), H(0), N(0) {}

Lattice::Lattice(const Settings &set, const Lattice::state_t &state):
    L(set.lattice.L), W(set.lattice.W), H(set.lattice.H)
{
    //SetFile("Lattice");
    //SetStream(&std::cout);
    N = L * W * H;
    periodic_L = set.lattice_boundary_conditions.periodic_boundary_condition_L;
    periodic_W = set.lattice_boundary_conditions.periodic_boundary_condition_W;
    periodic_H = set.lattice_boundary_conditions.periodic_boundary_condition_H;
    Construct();
    switch(state)
    {
        case Isotropic:
            IsotropicState();
            break;
        case IsotropicRighthanded:
            IsotropicRighthandedState();
            break;
        case Biaxial:
            BiaxialState();
            break;
        case BiaxialRighthanded:
            BiaxialRighthandedState();
            break;
        case BiaxialAlt:
            BiaxialStateAlt();
            break;
        case BiaxialRighthandedAlt:
            BiaxialRighthandedStateAlt();
            break;
    }
}

const Lattice &Lattice::operator=(const Lattice &s)
{
    L = s.L;
    W = s.W;
    H = s.H;
    N = L * W * H;
    periodic_L = s.periodic_L;
    periodic_W = s.periodic_W;
    periodic_H = s.periodic_H;
    Construct();
    for(int i = 0; i < N; i++)
        Particles[i].RestoreState(s.Particles[i]);
    return *this;
}

Lattice::Lattice(const Lattice &s)
{
    L = s.L;
    W = s.W;
    H = s.H;
    N = L * W * H;
    periodic_L = s.periodic_L;
    periodic_W = s.periodic_W;
    periodic_H = s.periodic_H;
    Construct();
    for(int i = 0; i < N; i++)
        Particles[i].RestoreState(s.Particles[i]);
}

Lattice::Lattice(Lattice &&s):
    L(s.L),
    W(s.W),
    H(s.H),
    N(s.N),
    periodic_L(s.periodic_L),
    periodic_W(s.periodic_W),
    periodic_H(s.periodic_H),
    Particles(std::move(s.Particles))
{

}

const Lattice &Lattice::operator=(Lattice && s)
{
    L = s.L;
    W = s.W;
    H = s.H;
    N = s.N;
    periodic_L = s.periodic_L;
    periodic_W = s.periodic_W;
    periodic_H = s.periodic_H;
    Particles = std::move(s.Particles);

    return *this;
}

void Lattice::IsotropicState()
{
    for(Particle & p : Particles)
    {
        p.SetOrientation(RandomPointOn4DSphereMarsaglia(1.0), plusminusone());
    }
}

void Lattice::IsotropicRighthandedState()
{
    for(Particle & p : Particles)
    {
        p.SetOrientation(RandomPointOn4DSphereMarsaglia(1.0), 1);
    }
}

void Lattice::BiaxialState()
{
    vect o(4);
    o[0] = 1.0;
    o[1] = o[2] = o[3] = 0.0;
    for(Particle & p : Particles)
    {
        p.SetOrientation(o, plusminusone());
    }
}

void Lattice::BiaxialRighthandedState()
{
    vect o(4);
    o[0] = 1.0;
    o[1] = o[2] = o[3] = 0.0;
    for(Particle & p : Particles)
    {
        p.SetOrientation(o, 1);
    }
}

void Lattice::BiaxialStateAlt()
{
    vect o(4);
    o[2] = o[3] = 0.0;
    o[0] = 1. / std::sqrt(2);
    o[1] = 1. / std::sqrt(2);
    for(Particle &p : Particles)
    {
        p.SetOrientation(o, plusminusone());
    }
}

void Lattice::BiaxialRighthandedStateAlt()
{
    vect o(4);
    o[2] = o[3] = 0.0;
    o[0] = 1. / std::sqrt(2);
    o[1] = 1. / std::sqrt(2);
    for(Particle &p : Particles)
    {
        p.SetOrientation(o, 1);
    }
}

void Lattice::Sweep(MCProto * proto, int &acc_rot, int &acc_p)
{
    for(int i = 0; i < N; i++)
    {
        int site = int(N * random01());
        int ar = 0, ap = 0;
        Particles[site].Nudge(proto, ar, ap);
        acc_rot += ar;
        acc_p += ap;

    }
}

const int &Lattice::GetN() const
{
    return N;
}

const int &Lattice::GetL() const
{
    return L;
}

const int &Lattice::GetW() const
{
    return W;
}

const int &Lattice::GetH() const
{
    return H;
}

const std::vector<Particle> &Lattice::GetParticles() const
{
    return Particles;
}

const double Lattice::GetMeanEPM() const
{
    double epm = 0.0;
    for(const Particle & p : Particles)
    {
        epm += p.GetEnergy();
    }
    return epm / N;
}

const bool &Lattice::GetPeriodicL() const
{
    return periodic_L;
}

const bool &Lattice::GetPeriodicW() const
{
    return periodic_W;
}

const bool &Lattice::GetPeriodicH() const
{
    return periodic_H;
}

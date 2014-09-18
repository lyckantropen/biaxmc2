/**
 * @file Lattice.h
 */

#ifndef _LATTICE_H
#define _LATTICE_H

#include "Particle.h"
#include "std.h"
#include "serializer.h"
#include "Settings.h"
#include "Random.h"
#include "4DSphereRW.h"
#include "ILoggable.h"

/**
 * @brief Stores the lattice
 *
 * This class stores the particles which constitute the lattice. It takes care
 * of creating and arranging the lattice and also performs the lattice-wide
 * Monte Carlo sweep.
 */

class Lattice: public ILoggable
{

    friend class PRE79StandardProperties;
    friend class SpatialCorrelation;
    friend class SpatialCorrelationEvolution;
    template<class stream_t>
    friend void operator|(boostbase::outserializer<stream_t> & s, Lattice & lat);
    template<class stream_t>
    friend void operator|(boostbase::inserializer<stream_t> & s, Lattice & lat);
    friend std::ostream & operator<<(std::ostream & s, const Lattice & lat);

    int N;      ///<number of particles
    int L;      ///<lattice length
    int W;      ///<lattice width
    int H;      ///<lattice height

    /// periodic boundary conditions flags
    bool periodic_L, periodic_W, periodic_H;

    std::vector<Particle>   Particles;          ///<linear vector of 'N' particles



    /// Construct the lattice. Invokes the Particle::Connect() function NxN times,
    /// while taking care of the periodic boundary conditions.
    ///
    /// @todo: Should the boundary conditions ever be changed, it needs to be done here.
    void Construct();
public:

    /// the many options for the initial conditions
    typedef enum { Isotropic, IsotropicRighthanded, Biaxial, BiaxialRighthanded, BiaxialAlt, BiaxialRighthandedAlt } state_t;

    /// Empty constructor for serialization purposes
    Lattice();


    /// Main constructor. Construct the lattice and initialize it with proper initial conditions.
    ///
    /// @param l Length
    /// @param w Width
    /// @param h Height
    /// @param state The desired option for the initial condition of the lattice
    Lattice(const Settings & set, const state_t & state = Isotropic);

    const Lattice & operator=(const Lattice & s);
    Lattice(const Lattice & s);
    Lattice(Lattice && s);
    const Lattice & operator=(Lattice && s);

private:
    /// Both degrees of freedom are randomly drawn
    void IsotropicState();

    /// The orientations are randomly drawn over the 4D-sphere, the parity is equal to +1
    void IsotropicRighthandedState();

    /// Ordered orientations, with the longest quadrupolar axis along Z, parity is random.
    void BiaxialState();

    /// Ordered orientations, with the longest quadrupolar axis along Z, parity is +1.
    void BiaxialRighthandedState();

    /// Ordered orientations, with the longest quadrupolar axis along X, parity is random.
    void BiaxialStateAlt();

    /// Ordered orientations, with the longest quadrupolar axis along X, parity is +1.
    void BiaxialRighthandedStateAlt();

public:
    /// A signle Monte Carlo sweep over all particles. The particles are picked at
    /// ranomd and updated using the Particle::Nudge() function. An implementation
    /// of the #MCProto prototype class is necessary.
    ///
    /// @param proto     An implementation of the #MCProto mechanism
    /// @param acc_rot   The counter which is incremented by +1 if a rotational move is accepted
    /// @param acc_p     The counter which is incremented by +1 if a parity move is accepted
    void Sweep(MCProto * proto, int & acc_rot, int & acc_p);

    ///Read-only accessor to #N
    const int & GetN() const;
    ///Read-only accessor to #L
    const int & GetL() const;
    ///Read-only accessor to #W
    const int & GetW() const;
    ///Read-only accessor to #H
    const int & GetH() const;
    ///Read-only accessor to #Particles
    const std::vector<Particle> & GetParticles() const;
    /**
     * Calculate the mean energy per molecule. The mean is the arithmetic mean.
     *
     * @return Arithmetic mean energy per molecule
     */
    const double GetMeanEPM() const;
};

template<class stream_t>
void operator|(boostbase::outserializer<stream_t> & s, Lattice & lat)
{
    //std::cout << "Saving\n";
    s | lat.N;
    s | lat.L;
    s | lat.W;
    s | lat.H;
    s | lat.periodic_L;
    s | lat.periodic_W;
    s | lat.periodic_H;
    for(const Particle & p : lat.GetParticles())
    {
        s | const_cast<Particle&>(p);
    }
}

template<class stream_t>
void operator|(boostbase::inserializer<stream_t> & s, Lattice & lat)
{
    //std::cout << "Loading\n";
    s | lat.N;
    s | lat.L;
    s | lat.W;
    s | lat.H;
    s | lat.periodic_L;
    s | lat.periodic_W;
    s | lat.periodic_H;
    lat.Construct();
    for(const Particle & p : lat.GetParticles())
    {
        s | const_cast<Particle&>(p);
    }
}

std::ostream & operator<<(std::ostream & s, const Lattice & lat);


#endif  /* _LATTICE_H */


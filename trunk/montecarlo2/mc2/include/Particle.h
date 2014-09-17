/**
 * @file Particle.h
 */

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "std.h"
#include "serializer.h"
#include "Maths.h"
#include "boost.h"
#include "MCProto.h"

#define coord_num 6
#define sqrt6   2.449489743

/**
 * @brief Holds information about the orientation of a single particle, also contains the atomic MC update
 *
 * This class holds the information about the orientation of a single particle
 * on the lattice. The information about the partilce's neighbors is stored.
 * Necessary quantities, such as the relevant molecular tensors are calculated
 * and stored.
 *
 * This class does not define interactions or the Monte Carlo method, which are
 * defined by implementing the #Hamiltonian and #MCProto prototype classes.
 */

class Particle
{
//    friend class SpatialCorrelationEvolution;                                   ///<ugly, should be rewritten
//    friend class Lattice;                                                       ///<ugly, should be rewritten
    ///serialization operator
    template<class serializer_t>
    friend void operator|(serializer_t & s, Particle & p);
    ///output operator
    friend std::ostream & operator<<(std::ostream & o, const Particle & p);
    short parity;       ///<parity degree of freedom, plus or minus one
    vect    x;          ///<quaternion which parametrizes the orientation of the molecule
    vect    ex;         ///<the X vector of the molecular orthonormal tripod
    vect    ey;         ///<the Y vector of the molecular orthonormal tripod
    vect    ez;         ///<the Z vector of the molecular orthonormal tripod
    vect    Qx;         ///<6 independent components of the X(x)X tensor
    vect    Qy;         ///<6 independent components of the Y(x)Y tensor
    vect    Qz;         ///<6 independent components of the Z(x)Z tensor
    vect    T;          ///<10 independent components of the T32 tensor
    vect    R;          ///<position on lattice
    double  energy;     ///<total energy per molecule
    std::vector<Particle*>  neighbors;  ///<6 pointers to the neighboring particles
    std::vector<int>        neighbors_indices;  ///<a list of 6 integers which specifiy their indices in Lattice::Particles

    /**
     * Update the energy of the particle using the two-particle energy between
     * all neighbors as well as the single-particle energy. The prototypes of the
     * 2 particle and single particle energy are declared in the #Hamiltonian
     * virtual function, but the actual parameter passed to this function needs
     * to link to an actual implementation of these functions.
     *
     * Thereby the interaction is NOT defined here.
     *
     * @param hamiltonian Pointer to the hamiltonian implementation.
     */
    void UpdateEnergy(const Hamiltonian * hamiltonian);
    /**
     * Update the energy of the neighboring particles. This necessarily implies
     * that a particle can have its energy updated many times during a sweep.
     * Again, to calculate the energy, we need the #Hamiltonian implementation
     * in the form of the parameter passed.
     *
     * @param h Pointer to the hamiltonian implementation
     */
    void UpdateNeighborsEnergy(const Hamiltonian * h);
public:
    ///Constructor
    Particle();
    /**
     * Copying constructor
     * @param s Source particle
     */
    Particle(const Particle & s);
    ///move constructor
    Particle(Particle && s);

    /**
     * Assignment operator
     * @param s Source particle
     * @return A copy of the assigned particle
     */
    const Particle & operator=(const Particle & s);

    ///move assignment operator
    const Particle & operator=(Particle && s);

    /**
     * Populate the lists #neighbors and #neighbors_indices by connecting the
     * specified particle to the present particle.
     *
     * @param n Particle to make connection with
     * @param i Index of the particle in the great particle container Lattice::Particles
     * @return False only if we attempt to connect more than coord_num particles
     */
    void Connect(Particle & n, const int & i);
    /**
     * A signle Monte Carlo update. Note that the actual acceptance scheme is not defined.
     * A prototype #MCProto exists, it however needs to be implemented further. An
     * example implementation is the #Metropolis class.
     *
     * The energies are updated initially and the resulting change of energy is
     * passed to the #MCProto implementation, which decides whether to accept or
     * reject the move. If the move is not accepted, the state is restored.
     *
     * The move is performed for the orientational and parity degrees
     * of freedom separately.
     *
     * The function also alters external counters for the rotational and parity
     * acceptance.
     *
     * @param proto
     * @param acc_rot   The parameter which is incremented by 1 if the rotational trial move is accepted
     * @param acc_p     The parameter which is incremented by 1 if the parity trial move is accepted
     */
    void Nudge(MCProto * proto, int & acc_rot, int & acc_p);

    /**
     * Restore all variables from another copy of a particle. Used primarily when a
     * MC move has been rejected and we need to go back.
     *
     * @param p A particle which we want to copy the variable values from
     *
     */
    void RestoreState(const Particle & p);

    /**
     * Set the orientation quaternion and the parity. Subsequently, all the
     * relevant tensors are calculated in this function
     *
     * @param X 4-component vector specifying the new orientation quaternion
     * @param p The new parity
     */
    void SetOrientation(const vect & X, const short & p);

    ///place at R = (x,y,z)
    void PlaceAt(int x, int y, int z);

    /**
     * Read only accessor to #energy
     * @return read only reference to #energy
     */
    const double & GetEnergy() const;
    /**
     * Read only accessor to #parity
     * @return read only reference to #parity
     */
    const short & GetParity() const;
    /**
     * Read only accessor to #x
     * @return read only reference to #x
     */
    const vect & GetX() const;
    /**
     * Read only accessor to #ex
     * @return read only reference to #ex
     */
    const vect & GetEX() const;
    /**
     * Read only accessor to #ey
     * @return read only reference to #ey
     */
    const vect & GetEY() const;
    /**
     * Read only accessor to #ez
     * @return read only reference to #ez
     */
    const vect & GetEZ() const;
    /**
     * Read only accessor to #Qx
     * @return read only reference to #Qx
     */
    const vect & GetQX() const;
    /**
     * Read only accessor to #Qy
     * @return read only reference to #Qy
     */
    const vect & GetQY() const;
    /**
     * Read only accessor to #Qz
     * @return read only reference to #Qz
     */
    const vect & GetQZ() const;
    /**
     * Read only accessor to #T
     * @return read only reference to #T
     */
    const vect & GetT() const;
    /**
     * Read only accessor to R
     * @return read only reference to R
     */
    const vect & GetR() const;

    int GetNNeighbors() const;
    int GetNeighborIndex(int k) const;
    const Particle *GetNeighbor(int k) const;

};
/**
 * Serializer operator
 * @param s
 * @param p
 */
template<class serializer_t>
void operator|(serializer_t & s, Particle & p)
{
    s | p.parity;
    s | p.x;
    s | p.ex;
    s | p.ey;
    s | p.ez;
    s | p.Qx;
    s | p.Qy;
    s | p.Qz;
    s | p.T;
    s | p.energy;
    s | p.neighbors_indices;
    s | p.R;
}
std::ostream & operator<<(std::ostream & o, const Particle & p);



#endif  /* _PARTICLE_H */


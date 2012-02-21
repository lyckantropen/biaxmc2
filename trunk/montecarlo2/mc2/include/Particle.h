/**
 * @file Particle.h
 */

#ifndef _PARTICLE_H
#define	_PARTICLE_H

#include "std.h"
#include "serializer.h"
#include "valarray_external.h"
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

class Particle {
    friend class SpatialCorrelationEvolution;                                   ///<ugly, should be rewritten
    friend class Lattice;                                                       ///<ugly, should be rewritten
    ///serialization operator
    template<class serializer_t>
    friend void operator|(serializer_t & s,Particle & p);
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
     * Restore all variables from another copy of a particle. Used primarily when a
     * MC move has been rejected and we need to go back.
     * 
     * @param p A particle which we want to copy the variable values from
     *
     */
    void RestoreState(const Particle & p){
        R=p.R;
        x=p.x;
        parity=p.parity;
        ex=p.ex;
        ey=p.ey;
        ez=p.ez;
        Qx=p.Qx;
        Qy=p.Qy;
        Qz=p.Qz;
        T=p.T;
        energy = p.energy;
    }
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
    void UpdateEnergy(shared_ptr<Hamiltonian> hamiltonian){
        if(hamiltonian==NULL) return;
        energy=0;
        for(int i=0;i<neighbors.size();i++){
            energy+=hamiltonian->TwoParticleEnergy(*this,*neighbors[i])/2.0;
        }
        energy+=hamiltonian->ExternalInteractionEnergy(*this);
    }
    /**
     * Update the energy of the neighboring particles. This necessarily implies
     * that a particle can have its energy updated many times during a sweep. 
     * Again, to calculate the energy, we need the #Hamiltonian implementation
     * in the form of the parameter passed.
     * 
     * @param h Pointer to the hamiltonian implementation
     */
    void UpdateNeighborsEnergy(shared_ptr<Hamiltonian> h){
        foreach(Particle * n, neighbors){
            n->UpdateEnergy(h);
        }
    }
    /**
     * Set the orientation quaternion and the parity. Subsequently, all the
     * relevant tensors are calculated in this function
     * 
     * @param X 4-component vector specifying the new orientation quaternion 
     * @param p The new parity
     */
    void SetOrientation(const vect & X, const short & p){
        parity = p;
        x=X;
        double  x11=x[1]*x[1],
                x22=x[2]*x[2],
                x33=x[3]*x[3],
                x01=x[0]*x[1],
                x02=x[0]*x[2],
                x03=x[0]*x[3],
                x12=x[1]*x[2],
                x13=x[1]*x[3],
                x23=x[2]*x[3];

        ex[0]=2*(-x22-x33+0.5);
	ex[1]=2*( x12+x03);
	ex[2]=2*( x13-x02);

	ey[0]=2*( x12-x03);
	ey[1]=2*(-x11-x33+0.5);
	ey[2]=2*( x01+x23);

	ez[0]=2*( x02+x13);
	ez[1]=2*(-x01+x23);
	ez[2]=2*(-x22-x11+0.5);

        Qx[0]=ex[0]*ex[0];  Qx[1]=ex[0]*ex[1];  Qx[2]=ex[0]*ex[2];
                            Qx[3]=ex[1]*ex[1];  Qx[4]=ex[1]*ex[2];
                                                Qx[5]=ex[2]*ex[2];

        Qy[0]=ey[0]*ey[0];  Qy[1]=ey[0]*ey[1];  Qy[2]=ey[0]*ey[2];
                            Qy[3]=ey[1]*ey[1];  Qy[4]=ey[1]*ey[2];
                                                Qy[5]=ey[2]*ey[2];

        Qz[0]=ez[0]*ez[0];  Qz[1]=ez[0]*ez[1];  Qz[2]=ez[0]*ez[2];
                            Qz[3]=ez[1]*ez[1];  Qz[4]=ez[1]*ez[2];
                                                Qz[5]=ez[2]*ez[2];

        //000
        T[0]=6.0*  (ex[0]*ey[0]*ez[0]);     // 1
        //100
        T[1]=2.0*  (ex[0]*ey[0]*ez[1]+      // 3
                    ex[0]*ey[1]*ez[0]+
                    ex[1]*ey[0]*ez[0]);
        //110
        T[2]=2.0*  (ex[0]*ey[1]*ez[1]+      // 3
                    ex[1]*ey[0]*ez[1]+
                    ex[1]*ey[1]*ez[0]);
        //111
        T[3]=6.0*  (ex[1]*ey[1]*ez[1]);     // 1
        //200
        T[4]=2.0*  (ex[0]*ey[0]*ez[2]+      // 3
                    ex[0]*ey[2]*ez[0]+
                    ex[2]*ey[0]*ez[0]);
        //210
        T[5]=      (ex[0]*ey[1]*ez[2]+      // 6
                    ex[0]*ey[2]*ez[1]+
                    ex[1]*ey[0]*ez[2]+
                    ex[2]*ey[0]*ez[1]+
                    ex[1]*ey[2]*ez[0]+
                    ex[2]*ey[1]*ez[0]
                    );
        //211
        T[6]=2.0*  (ex[1]*ey[1]*ez[2]+      // 3
                    ex[1]*ey[2]*ez[1]+
                    ex[2]*ey[1]*ez[1]
                    );
        //220
        T[7]=2.0*  (ex[2]*ey[2]*ez[0]+      // 3
                    ex[0]*ey[2]*ez[2]+
                    ex[2]*ey[0]*ez[2]);
        //221
        T[8]=2.0*  (ex[2]*ey[2]*ez[1]+      // 3
                    ex[1]*ey[2]*ez[2]+
                    ex[2]*ey[1]*ez[2]);
        //222
        T[9]=6.0*  (ex[2]*ey[2]*ez[2]);     // 1

        T*=(parity/sqrt6);
        
        //T/=sqrt6;
    }
public:
    ///Constructor
    Particle():
    x(4),ex(3),ey(3),ez(3),Qx(6),Qy(6),Qz(6),T(10),R(3) {
        //,e2p(coord_num){
        parity=1;
        energy=0;
    }
    /**
     * Copying constructor
     * @param s Source particle
     */
    Particle(const Particle & s){
        x.resize(4,0.0);
        ex.resize(3,0.0);
        ey.resize(3,0.0);
        ez.resize(3,0.0);
        Qx.resize(6,0.0);
        Qy.resize(6,0.0);
        Qz.resize(6,0.0);
        T.resize(10,0.0);
        R.resize(3,0.0);
        x=s.x;
        ex=s.ex;
        ey=s.ey;
        ez=s.ez;
        Qx=s.Qx;
        Qy=s.Qy;
        Qz=s.Qz;
        T=s.T;
        R=s.R;
        parity=s.parity;
        energy=s.energy;
        neighbors = s.neighbors;
        neighbors_indices = s.neighbors_indices;        
    }
    /**
     * Assignment operator
     * @param s Source particle
     * @return A copy of the assigned particle
     */
    const Particle & operator=(const Particle & s){
        x.resize(4,0.0);
        ex.resize(3,0.0);
        ey.resize(3,0.0);
        ez.resize(3,0.0);
        Qx.resize(6,0.0);
        Qy.resize(6,0.0);
        Qz.resize(6,0.0);
        T.resize(10,0.0);
        R.resize(3,0.0);
        R=s.R;
        x=s.x;
        ex=s.ex;
        ey=s.ey;
        ez=s.ez;
        Qx=s.Qx;
        Qy=s.Qy;
        Qz=s.Qz;
        T=s.T;
        parity=s.parity;
        energy=s.energy;
        neighbors = s.neighbors;
        neighbors_indices = s.neighbors_indices;        
        return * this;
    }
    /**
     * Populate the lists #neighbors and #neighbors_indices by connecting the
     * specified particle to the present particle.
     * 
     * @param n Particle to make connection with
     * @param i Index of the particle in the great particle container Lattice::Particles
     * @return False only if we attempt to connect more than coord_num particles
     */
    bool Connect(Particle & n, const int & i){
        if(neighbors.size()<coord_num){
            neighbors.push_back(&n);
            neighbors_indices.push_back(i);
            return true;
        }
        else return false;
    }
    /**
     * A signle Monte Carlo update. Note that the actual acceptance scheme is not defined.
     * A prototype #MCProto exists, it however needs to be implemented further. An
     * example implementation is the #Metropolis class.
     * 
     * The energies are updated initially and the resulting change of energy is
     * passed to the #MCProto implementation, which decides whether to accept or
     * reject the move. If the move is not accepted, the state is restored.
     * 
     * The move is performed separately for the orientational and parity degrees
     * of freedom separately.
     * 
     * The function also alters external counters for the rotational and parity
     * acceptance.
     * 
     * @param proto
     * @param acc_rot   The parameter which is incremented by 1 if the rotational trial move is accepted
     * @param acc_p     The parameter which is incremented by 1 if the parity trial move is accepted
     */
    void Nudge(shared_ptr<MCProto> proto, int & acc_rot, int & acc_p){
        Particle old_state = *this;
        vect newX = proto->OrientationNudge(x);
        short newP = proto->ParityNudge(parity);


        //SetOrientation(proto->OrientationNudge(x),proto->ParityNudge(parity));
        
        //ruch rotacyjny
        SetOrientation(newX,parity);
        UpdateEnergy(proto->GetHamiltonian());
        // mnożymy przez 2 bo mamy energię na cząstkę
        if(!proto->Accept(2.0*(energy-old_state.energy))){
            // ruch niezaakceptowany
            RestoreState(old_state);
            acc_rot+=0;
        }
        else{
            // ruch zaakceptowany
            UpdateNeighborsEnergy(proto->GetHamiltonian());
            acc_rot+=1;
        }
        //--

        old_state = *this;
        //ruch parzystości
        SetOrientation(x,newP);
        UpdateEnergy(proto->GetHamiltonian());
        // mnożymy przez 2 bo mamy energię na cząstkę
        if(!proto->Accept(2.0*(energy-old_state.energy))){
            // ruch niezaakceptowany
            RestoreState(old_state);
            acc_p+=0;
        }
        else{
            // ruch zaakceptowany
            UpdateNeighborsEnergy(proto->GetHamiltonian());
            acc_p+=1;
        }
        //--
    }

    /**
     * Read only accessor to #energy
     * @return read only reference to #energy
     */
    const double & GetEnergy() const {
        return energy;
    }
    /**
     * Read only accessor to #parity
     * @return read only reference to #parity
     */
    const short & GetParity() const {
        return parity;
    }
    /**
     * Read only accessor to #x
     * @return read only reference to #x
     */
    const vect & GetX() const {
        return x;
    }
    /**
     * Read only accessor to #ex
     * @return read only reference to #ex
     */
    const vect & GetEX() const {
        return ex;
    }
    /**
     * Read only accessor to #ey
     * @return read only reference to #ey
     */
    const vect & GetEY() const {
        return ey;
    }
    /**
     * Read only accessor to #ez
     * @return read only reference to #ez
     */
    const vect & GetEZ() const {
        return ez;
    }
    /**
     * Read only accessor to #Qx
     * @return read only reference to #Qx
     */
    const vect & GetQX() const {
        return Qx;
    }
    /**
     * Read only accessor to #Qy
     * @return read only reference to #Qy
     */
    const vect & GetQY() const {
        return Qy;
    }
    /**
     * Read only accessor to #Qz
     * @return read only reference to #Qz
     */
    const vect & GetQZ() const {
        return Qz;
    }
    /**
     * Read only accessor to #T
     * @return read only reference to #T
     */
    const vect & GetT() const {
        return T;
    }
    /**
     * Read only accessor to R
     * @return read only reference to R
     */
    const vect & GetR() const {
        return R;
    }

};
/**
 * Serializer operator
 * @param s
 * @param p
 */
template<class serializer_t>
void operator|(serializer_t & s,Particle & p){
    s|p.parity;
    s|p.x;
    s|p.ex;
    s|p.ey;
    s|p.ez;
    s|p.Qx;
    s|p.Qy;
    s|p.Qz;
    s|p.T;
    s|p.energy;
    s|p.neighbors_indices;
    s|p.R;
}
extern std::ostream & operator<<(std::ostream & o, const Particle & p);



#endif	/* _PARTICLE_H */


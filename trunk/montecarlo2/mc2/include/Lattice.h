/**
 * @file Lattice.h
 */

#ifndef _LATTICE_H
#define	_LATTICE_H

#include "Particle.h"
#include "std.h"
#include "serializer.h"
#include "Settings.h"
#include "random01.h"

/**
 * @brief Stores the lattice
 * 
 * This class stores the particles which constitute the lattice. It takes care 
 * of creating and arranging the lattice and also performs the lattice-wide 
 * Monte Carlo sweep.
 */

class Lattice {
    friend class PRE79StandardProperties;
    friend class SpatialCorrelation;
    friend class SpatialCorrelationEvolution;
    template<class stream_t>
    friend void operator|(boostbase::outserializer<stream_t> & s, Lattice & lat);
    template<class stream_t>
    friend void operator|(boostbase::inserializer<stream_t> & s, Lattice & lat);
    friend std::ostream & operator<<(std::ostream & s,const Lattice & lat);
    int N;      ///<number of particles
    int L;      ///<lattice length
    int W;      ///<lattice width
    int H;      ///<lattice height
    bool periodic_L, periodic_W, periodic_H;
    std::vector<Particle>   Particles;          ///<linear vector of #N particles
    /**
     * Construct the lattice. Invokes the Particle::Connect() function #N*12 times,
     * while taking care of the periodic boundary conditions.
     * 
     * @todo: Should the boundary conditions ever be changed, it needs to be done here.
     */
    void Construct(){
        //std::cout << "Constructing\n";
        if(Particles.size()!=0) Particles.clear();
        Particles.resize(N); //= std::vector<Particle>(N);
        
        for(int i=0;i<N;i++){
            Particle & cur = Particles[i];
            int A=L*W;
            int x=0,y=0,z=0;
            x=i%L;
            y=i%A/L;
            z=i/A;
            cur.R[0]=x;
            cur.R[1]=y;
            cur.R[2]=z;
            
            if(x==L-1 && periodic_L)
                cur.Connect(Particles[i-L+1],i-L+1);
            else
                if(L!=1) cur.Connect(Particles[i+1],i+1);
            if(x==0 && periodic_L)
                cur.Connect(Particles[i+L-1],i+L-1);
            else
                if(L!=1) cur.Connect(Particles[i-1],i-1);
            
            if(y==W-1 && periodic_W)
                cur.Connect(Particles[i-A+L],i-A+L);
            else
                if(W!=1) cur.Connect(Particles[i+L],i+L);
            if(y==0 && periodic_W)
                cur.Connect(Particles[i+A-L],i+A-L);
            else
                if(W!=1) cur.Connect(Particles[i-L],i-L);
            
            if(z==H-1 && periodic_H)
                cur.Connect(Particles[i-N+A],i-N+A);
            else
                if(H!=1) cur.Connect(Particles[i+A],i+A);
            if(z==0 && periodic_H)
                cur.Connect(Particles[i+N-A],i+N-A);
            else
                if(H!=1) cur.Connect(Particles[i-A],i-A);
        }

    }
public:
    ///the multiple options for the initial conditions
    typedef enum { Isotropic, IsotropicRighthanded, Biaxial, BiaxialRighthanded, BiaxialAlt, BiaxialRighthandedAlt } state_t;
    Lattice():L(0),W(0),H(0),N(0) {}
    /**
     * Main constructor. Construct the lattice and initialize it with proper initial conditions.
     * 
     * @param l Length
     * @param w Width
     * @param h Height
     * @param state The desired option for the initial condition of the lattice
     */
    Lattice(const Settings & set,const state_t & state=Isotropic):
    L(set.lattice.L),W(set.lattice.W),H(set.lattice.H)
    {
        N=L*W*H;
        periodic_L=set.lattice_boundary_conditions.periodic_boundary_condition_L;
        periodic_W=set.lattice_boundary_conditions.periodic_boundary_condition_W;
        periodic_H=set.lattice_boundary_conditions.periodic_boundary_condition_H;
        Construct();
        switch(state){
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

    /**
     * Assignment operator
     * 
     * @param s source lattice
     * @return 
     */
    const Lattice & operator=(const Lattice & s){
        L=s.L;
        W=s.W;
        H=s.H;
        N=L*W*H;
        periodic_L=s.periodic_L;
        periodic_W=s.periodic_W;
        periodic_H=s.periodic_H;
        Construct();
        for(int i=0;i<N;i++)
            Particles[i].RestoreState(s.Particles[i]);
        return *this;
    }
    /**
     * Copy constructor
     * @param s source lattice
     */
    Lattice(const Lattice & s){
        L=s.L;
        W=s.W;
        H=s.H;
        N=L*W*H;
        periodic_L=s.periodic_L;
        periodic_W=s.periodic_W;
        periodic_H=s.periodic_H;
        Construct();
        for(int i=0;i<N;i++)
            Particles[i].RestoreState(s.Particles[i]);
    }

private:
    /**
     * Both degrees of freedom are randomly drawn
     */
    void IsotropicState(){
        foreach(Particle & p, Particles){
            p.SetOrientation(RandomPointOn4DSphereMarsaglia(1.0),plusminusone());
        }
    }
    /**
     * The orientations are randomly drawn over the 4D-sphere, the parity is equal to +1
     */
    void IsotropicRighthandedState(){
        foreach(Particle & p, Particles){
            p.SetOrientation(RandomPointOn4DSphereMarsaglia(1.0),1);
        }
    }
    /**
     * Ordered orientations, with the longest quadrupolar axis along Z, parity is random.
     */
    void BiaxialState(){
        vect o(4);
        o[0]=1.0;
        o[1]=o[2]=o[3]=0.0;
        foreach(Particle & p, Particles){
            p.SetOrientation(o,plusminusone());
        }
    }
    /**
     * Ordered orientations, with the longest quadrupolar axis along Z, parity is +1.
     */
    void BiaxialRighthandedState(){
        vect o(4);
        o[0]=1.0;
        o[1]=o[2]=o[3]=0.0;
        foreach(Particle & p, Particles){
            p.SetOrientation(o,1);
        }
    }
        /**
     * Ordered orientations, with the longest quadrupolar axis along X, parity is random.
     */
    void BiaxialStateAlt(){
	    vect o(4);
	    o[2]=o[3]=0.0;
	    o[0]=1./std::sqrt(2);
	    o[1]=1./std::sqrt(2);
	    foreach(Particle &p, Particles){
		    p.SetOrientation(o,plusminusone());
	    }
    }
        /**
     * Ordered orientations, with the longest quadrupolar axis along X, parity is +1.
     */
    void BiaxialRighthandedStateAlt(){
	    vect o(4);
	    o[2]=o[3]=0.0;
	    o[0]=1./std::sqrt(2);
	    o[1]=1./std::sqrt(2);
	    foreach(Particle &p, Particles){
		    p.SetOrientation(o,1);
	    }
    }

public:
    /**
     * A signle Monte Carlo sweep over all particles. The particles are picked at
     * ranomd and updated using the Particle::Nudge() function. An implementation
     * of the #MCProto prototype class is necessary.
     * 
     * @param proto     An implementation of the #MCProto mechanism
     * @param acc_rot   The counter which is incremented by +1 if a rotational move is accepted
     * @param acc_p     The counter which is incremented by +1 if a parity move is accepted
     */
    void Sweep(shared_ptr<MCProto> proto, int & acc_rot, int & acc_p){
        for(int i=0;i<N;i++){
            int site = int(N*random01());
            int ar=0,ap=0;
            Particles[site].Nudge(proto,ar,ap);
            acc_rot+=ar;
            acc_p+=ap;
                
        }
    }

    ///Read-only accessor to #N
    const int & GetN() const {
        return N;
    }
    ///Read-only accessor to #L
    const int & GetL() const {
        return L;
    }
    ///Read-only accessor to #W
    const int & GetW() const {
        return W;
    }
    ///Read-only accessor to #H
    const int & GetH() const {
        return H;
    }
    ///Read-only accessor to #Particles
    const std::vector<Particle> & GetParticles() const {
        return Particles;
    }
    /**
     * Calculate the mean energy per molecule. The mean is the arithmetic mean.
     * 
     * @return Arithmetic mean energy per molecule 
     */
    const double GetMeanEPM() const {
        double epm=0.0;
        foreach(const Particle & p, Particles) {
            epm+=p.GetEnergy();
        }
        return epm/N;
    }
};

template<class stream_t>
void operator|(boostbase::outserializer<stream_t> & s, Lattice & lat){
    //std::cout << "Saving\n";
    s|lat.N;
    s|lat.L;s|lat.W;s|lat.H;
    s|lat.periodic_L;
    s|lat.periodic_W;
    s|lat.periodic_H;
    foreach(Particle & p, lat.Particles){
        s|p;
    }
}

template<class stream_t>
void operator|(boostbase::inserializer<stream_t> & s, Lattice & lat){
    //std::cout << "Loading\n";
    s|lat.N;
    s|lat.L;s|lat.W;s|lat.H;
    s|lat.periodic_L;
    s|lat.periodic_W;
    s|lat.periodic_H;
    lat.Construct();
    foreach(Particle & p, lat.Particles){
        s|p;
    }
}

extern std::ostream & operator<<(std::ostream & s,const Lattice & lat);


#endif	/* _LATTICE_H */


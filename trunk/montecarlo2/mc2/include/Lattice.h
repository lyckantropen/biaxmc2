/* 
 * File:   Lattice.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 17:41
 */

#ifndef _LATTICE_H
#define	_LATTICE_H

#include "Particle.h"
#include "std.h"
#include "serializer.h"
#include "random01.h"

class Lattice {
    friend class PRE79StandardProperties;
    friend class SpatialCorrelation;
    friend class SpatialCorrelationEvolution;
    template<class stream_t>
    friend void operator|(boostbase::outserializer<stream_t> & s, Lattice & lat);
    template<class stream_t>
    friend void operator|(boostbase::inserializer<stream_t> & s, Lattice & lat);
    friend std::ostream & operator<<(std::ostream & s,const Lattice & lat);
    int N;      ///<liczba cząstek
    int L,W,H;  ///<wymiary siatki
    std::vector<Particle>   Particles;
    void Construct(){
        //std::cout << "Constructing\n";
        if(Particles.size()!=0) Particles.clear();
        Particles = std::vector<Particle>(N);
        
        for(int i=0;i<N;i++){
            Particle & cur = Particles[i];
            int A=L*W;
            int x=0,y=0,z=0;
            x=i%L;
            y=i%A/L;
            z=i/A;
            if(x==L-1)
                cur.Connect(Particles[i-L+1],i-L+1);
            else
                cur.Connect(Particles[i+1],i+1);
            if(x==0)
                cur.Connect(Particles[i+L-1],i+L-1);
            else
                cur.Connect(Particles[i-1],i-1);
            if(y==W-1)
                cur.Connect(Particles[i-A+L],i-A+L);
            else
                cur.Connect(Particles[i+L],i+L);
            if(y==0)
                cur.Connect(Particles[i+A-L],i+A-L);
            else
                cur.Connect(Particles[i-L],i-L);
            if(z==H-1)
                cur.Connect(Particles[i-N+A],i-N+A);
            else
                cur.Connect(Particles[i+A],i+A);
            if(z==0)
                cur.Connect(Particles[i+N-A],i+N-A);
            else
                cur.Connect(Particles[i-A],i-A);
        }

    }
public:
    typedef enum { Isotropic, IsotropicRighthanded, Biaxial, BiaxialRighthanded } state_t;
    Lattice():L(0),W(0),H(0),N(0) {}
    Lattice(const int & l, const int & w, const int & h,const state_t & state=Isotropic):
    L(l),W(w),H(h),N(l*w*h)
    {
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
        }
    }

    // tak to trzeba zrobić, bo Particle zawiera wskaźniki
    const Lattice & operator=(const Lattice & s){
        L=s.L;
        W=s.W;
        H=s.H;
        N=L*W*H;
        Construct();
        for(int i=0;i<N;i++)
            Particles[i].RestoreState(s.Particles[i]);
        return *this;
    }
    Lattice(const Lattice & s){
        L=s.L;
        W=s.W;
        H=s.H;
        N=L*W*H;
        Construct();
        for(int i=0;i<N;i++)
            Particles[i].RestoreState(s.Particles[i]);
    }

private:
    ///najwyższa symetria
    void IsotropicState(){
        foreach(Particle & p, Particles){
            p.SetOrientation(RandomPointOn4DSphereMarsaglia(1.0),plusminusone());
        }
    }
    ///porządek chiralny
    void IsotropicRighthandedState(){
        foreach(Particle & p, Particles){
            p.SetOrientation(RandomPointOn4DSphereMarsaglia(1.0),1);
        }
    }
    ///porządek dwuosiowy
    void BiaxialState(){
        vect o(4);
        o[0]=1.0;
        o[1]=o[2]=o[3]=0.0;
        foreach(Particle & p, Particles){
            p.SetOrientation(o,plusminusone());
        }
    }
    ///porządek dwuosiowy + chiralny
    void BiaxialRighthandedState(){
        vect o(4);
        o[0]=1.0;
        o[1]=o[2]=o[3]=0.0;
        foreach(Particle & p, Particles){
            p.SetOrientation(o,1);
        }
    }
public:
    ///jedno przemiecenie Monte Carlo
    int Sweep(MCProto * proto){
        int accepted = 0;
        for(int i=0;i<N;i++){
            int site = int(N*random01());
            if(Particles[i].Nudge(proto))
                accepted++;
        }
        return accepted;
    }

    //Accessors
    const int & GetN() const {
        return N;
    }
    const int & GetL() const {
        return L;
    }
    const int & GetW() const {
        return W;
    }
    const int & GetH() const {
        return H;
    }
    const std::vector<Particle> & GetParticles() const {
        return Particles;
    }
};

template<class stream_t>
void operator|(boostbase::outserializer<stream_t> & s, Lattice & lat){
    //std::cout << "Saving\n";
    s|lat.N;
    s|lat.L;s|lat.W;s|lat.H;
    foreach(Particle & p, lat.Particles){
        s|p;
    }
}

template<class stream_t>
void operator|(boostbase::inserializer<stream_t> & s, Lattice & lat){
    //std::cout << "Loading\n";
    s|lat.N;
    s|lat.L;s|lat.W;s|lat.H;
    lat.Construct();
    foreach(Particle & p, lat.Particles){
        s|p;
    }
}

extern std::ostream & operator<<(std::ostream & s,const Lattice & lat);


#endif	/* _LATTICE_H */


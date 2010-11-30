/* 
 * File:   particle.h
 * Author: karol
 *
 * Created on 16 listopad 2009, 17:28
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

class Particle {
    friend class SpatialCorrelationEvolution;
    friend class Lattice;
    template<class serializer_t>
    friend void operator|(serializer_t & s,Particle & p);
    friend std::ostream & operator<<(std::ostream & o, const Particle & p);
    short parity;
    vect    x;          ///<orientacyjne stopnie swobody
    vect    ex,ey,ez;   ///<baza związana z molekułą
    vect    Qx,Qy,Qz;   ///<6 składników tensorów ei(x)ei, z nich można skonstruować inne tensory
    vect    T;          ///<10 składników tensora T32
    double  energy;     ///<całkowita energia 1-cząstkowa
    std::vector<Particle*>  neighbors;  ///<wskaźniki do sąsiadów
    std::vector<int>        neighbors_indices;  ///<indeksy sąsiadów (pomocnicze)

    void RestoreState(const Particle & p){
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
    void UpdateEnergy(Hamiltonian * hamiltonian){
        if(hamiltonian==NULL) return;
        energy=0;
        for(int i=0;i<neighbors.size();i++){
            energy+=hamiltonian->TwoParticleEnergy(*this,*neighbors[i])/2.0;
        }
        energy+=hamiltonian->ExternalInteractionEnergy(*this);
    }
    void UpdateNeighborsEnergy(Hamiltonian * h){
        foreach(Particle * n, neighbors){
            n->UpdateEnergy(h);
        }
    }
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
    Particle():
    x(4),ex(3),ey(3),ez(3),Qx(6),Qy(6),Qz(6),T(10) {
        //,e2p(coord_num){
        parity=1;
        energy=0;
    }
    bool Connect(Particle & n, const int & i){
        if(neighbors.size()<coord_num){
            neighbors.push_back(&n);
            neighbors_indices.push_back(i);
            return true;
        }
        else return false;
    }
    ///pojedynczy ruch Monte Carlo względem danego prototypu procesu i hamiltonianu
    void Nudge(MCProto * proto, int & acc_rot, int & acc_p){
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

    //Accessors
    const double & GetEnergy() const {
        return energy;
    }
    const short & GetParity() const {
        return parity;
    }
    const vect & GetX() const {
        return x;
    }
    const vect & GetEX() const {
        return ex;
    }
    const vect & GetEY() const {
        return ey;
    }
    const vect & GetEZ() const {
        return ez;
    }
    const vect & GetQX() const {
        return Qx;
    }
    const vect & GetQY() const {
        return Qy;
    }
    const vect & GetQZ() const {
        return Qz;
    }
    const vect & GetT() const {
        return T;
    }

};
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
}
extern std::ostream & operator<<(std::ostream & o, const Particle & p);



#endif	/* _PARTICLE_H */


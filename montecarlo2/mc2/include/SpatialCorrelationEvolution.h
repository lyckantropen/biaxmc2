/* 
 * File:   SpatialCorrelationEvolution.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 11:37
 */

#ifndef _SPATIALCORRELATIONEVOLUTION_H
#define	_SPATIALCORRELATIONEVOLUTION_H

//#include "SpatialCorrelation.h"
#include "Lattice.h"
#include "Statistical.h"
#include "Contractions.h"
#include "serializer.h"
#include "std.h"

/**
 * klasa prototypowa do obliczania funkcji korelacji
 */

class SpatialCorrelationEvolution {
public:
    typedef std::valarray<std::valarray<std::valarray<double> > >  vectijk;
private:
    template<class serializer_t>
    friend void operator|(serializer_t &,SpatialCorrelationEvolution &);
    std::vector<vect> correlation;    ///<ewolucja funkcji korelacji
    Lattice *   lat;
    int ncycles;
    int acc_idx;
    int max;
    bool readonly;
    vectijk nnn;    ///nnn[site][neighbor][distance], tablica indeksów sąsiadów
    void ConstructNNN(Lattice * lat){
        if(readonly) return;
        for(int site=0;site<nnn.size();site++)
            for(int n=0;n<coord_num;n++)
                for(int dist=0;dist<(max+1);dist++){
                    if(dist==0)
                        nnn[site][n][dist]=site;
                    else if(dist==1)
                        nnn[site][n][dist]=lat->Particles[site].neighbors_indices[n];
                    else
                        nnn[site][n][dist]=lat->Particles[nnn[site][n][dist-1]].neighbors_indices[n];
                }
    }

protected:
    virtual double CalculateContraction(const Particle &,const Particle&)   { return 0.0; };
private:
    virtual void CalculateCorrelation(){
        vect result(max+1);
        for(int site=0;site<lat->GetN();site++)
            for(int n=0;n<coord_num;n++)
                for(int dist=0;dist<(max+1);dist++){
                    result[dist]+=CalculateContraction(lat->Particles[site],lat->Particles[nnn[site][n][dist]]);
                }
        correlation[acc_idx]=result/double(lat->GetN()*coord_num);
    }
protected:
    ///konstruktor tradycyjny (zasłonięty, bo funkcja CalculateContraction pozostaje do określenia w klasie dziedziczącej)
    SpatialCorrelationEvolution(Lattice * l=NULL, int nc=0){
        if(l==NULL) return; //serializacja

        readonly=false;
        lat=l;
        max=lat->GetL()/2;
        ncycles=nc;
        acc_idx=-1;
        for(int i=0;i<ncycles;i++)
            correlation.push_back(vect(max+1));

        nnn.resize(lat->GetN());
        for(int i=0;i<nnn.size();i++){
            nnn[i].resize(6);
            for(int j=0;j<6;j++)
                nnn[i][j].resize(max+1);
        }
        ConstructNNN(lat);
    }
    template<class serializer_t>
    void SerializerHelper(serializer_t & s){
        s|ncycles;
        s|acc_idx;
        s|max;
        s|correlation;
    }
public:
    ///konstruktor serializacyjny
    SpatialCorrelationEvolution(){
        lat=NULL;
        ncycles=0;
        acc_idx=0;
        readonly=true;
    }
    ///policzenie kolejnej wartości funkcji korelacji
    void Update(){
        if(readonly) return;
        if((acc_idx+1)>=ncycles) return;
        acc_idx++;
        CalculateCorrelation();
    }
    ///granica funkcji korelacji dla maksymalnej odległości (do par. porządku)
    Value   Limit() const {
        vect cor(acc_idx+1);
        for(int i=0;i<(acc_idx+1);i++){
            cor[i]=correlation[i][max];
        }
        return BootstrapMean(cor);
    }
    
    ///dołożenie danych z innej ewolucji tej samej funkcji korelacji, nie należy wykonywać przed zakończeniem obliczania bieżącej funkcji
    void Append(const SpatialCorrelationEvolution & c){
        ncycles+=c.ncycles;
        acc_idx+=c.ncycles;
        correlation.insert(correlation.end(),c.correlation.begin(),c.correlation.end());
    }
    ///średnia z funkcji korelacji po czasie
    vect    Mean(int acc_idx=0) const {
        return MeanVector(correlation,0,acc_idx);
    }
    //Accessors
    const int & GetNCycles() const {
        return ncycles;
    }
    const int & GetAccIdx() const {
        return acc_idx;
    }
    const vectijk & GetNNN() const {
        return nnn;
    }
    const vect & CurrentCorrelation() const {
        return correlation[acc_idx];
    }
    const int & GetMax() const {
        return max;
    }
    const vect & operator[](int t) const {
        return correlation[t];
    }
};

template<class serializer_t>
void operator|(serializer_t & s,SpatialCorrelationEvolution & e){
    e.SerializerHelper<serializer_t>(s);
};

#endif	/* _SPATIALCORRELATIONEVOLUTION_H */


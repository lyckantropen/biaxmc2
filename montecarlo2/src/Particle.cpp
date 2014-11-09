#include "Particle.h"

std::ostream & operator<<(std::ostream & s, const Particle & p)
{
    s << p.parity ;
    for(int i = 0; i < p.x.size(); i++)
    {
        s << " " << p.x[i] ;
    }
    return s;
}



void Particle::PlaceAt(int x, int y, int z)
{
    R = { static_cast<double>(x) , static_cast<double>(y) , static_cast<double>(z) };
}

const double &Particle::GetEnergy() const
{
    return energy;
}

const short &Particle::GetParity() const
{
    return parity;
}

const vect &Particle::GetX() const
{
    return x;
}

const vect &Particle::GetEX() const
{
    return ex;
}

const vect &Particle::GetEY() const
{
    return ey;
}

const vect &Particle::GetEZ() const
{
    return ez;
}

const vect &Particle::GetQX() const
{
    return Qx;
}

const vect &Particle::GetQY() const
{
    return Qy;
}

const vect &Particle::GetQZ() const
{
    return Qz;
}

const vect &Particle::GetT() const
{
    return T;
}

const vect &Particle::GetR() const
{
    return R;
}

int Particle::GetNNeighbors() const
{
    return neighbors.size();
}

int Particle::GetNeighborIndex(int k) const
{
    return neighbors_indices[k];
}

const Particle * Particle::GetNeighbor(int k) const
{
    return neighbors[k];
}

void Particle::Connect(Particle &n, const int &i)
{
    //std::cout << "Connecting to " << i << "\n";
    neighbors.push_back(&n);
    neighbors_indices.push_back(i);

}


void Particle::RestoreState(const Particle &p)
{
    R = p.R;
    x = p.x;
    parity = p.parity;
    ex = p.ex;
    ey = p.ey;
    ez = p.ez;
    Qx = p.Qx;
    Qy = p.Qy;
    Qz = p.Qz;
    T = p.T;
    energy = p.energy;
}





void Particle::SetOrientation(const vect &X, const short &p)
{
    parity = p;
    x = X;
    double  x11 = x[1] * x[1],
            x22 = x[2] * x[2],
            x33 = x[3] * x[3],
            x01 = x[0] * x[1],
            x02 = x[0] * x[2],
            x03 = x[0] * x[3],
            x12 = x[1] * x[2],
            x13 = x[1] * x[3],
            x23 = x[2] * x[3];

    ex[0] = 2 * (-x22 - x33 + 0.5);
    ex[1] = 2 * (x12 + x03);
    ex[2] = 2 * (x13 - x02);

    ey[0] = 2 * (x12 - x03);
    ey[1] = 2 * (-x11 - x33 + 0.5);
    ey[2] = 2 * (x01 + x23);

    ez[0] = 2 * (x02 + x13);
    ez[1] = 2 * (-x01 + x23);
    ez[2] = 2 * (-x22 - x11 + 0.5);

    Qx[0] = ex[0] * ex[0];
    Qx[1] = ex[0] * ex[1];
    Qx[2] = ex[0] * ex[2];
    Qx[3] = ex[1] * ex[1];
    Qx[4] = ex[1] * ex[2];
    Qx[5] = ex[2] * ex[2];

    Qy[0] = ey[0] * ey[0];
    Qy[1] = ey[0] * ey[1];
    Qy[2] = ey[0] * ey[2];
    Qy[3] = ey[1] * ey[1];
    Qy[4] = ey[1] * ey[2];
    Qy[5] = ey[2] * ey[2];

    Qz[0] = ez[0] * ez[0];
    Qz[1] = ez[0] * ez[1];
    Qz[2] = ez[0] * ez[2];
    Qz[3] = ez[1] * ez[1];
    Qz[4] = ez[1] * ez[2];
    Qz[5] = ez[2] * ez[2];

    //000
    T[0] = 6.0 * (ex[0] * ey[0] * ez[0]); // 1
    //100
    T[1] = 2.0 * (ex[0] * ey[0] * ez[1] + // 3
                  ex[0] * ey[1] * ez[0] +
                  ex[1] * ey[0] * ez[0]);
    //110
    T[2] = 2.0 * (ex[0] * ey[1] * ez[1] + // 3
                  ex[1] * ey[0] * ez[1] +
                  ex[1] * ey[1] * ez[0]);
    //111
    T[3] = 6.0 * (ex[1] * ey[1] * ez[1]); // 1
    //200
    T[4] = 2.0 * (ex[0] * ey[0] * ez[2] + // 3
                  ex[0] * ey[2] * ez[0] +
                  ex[2] * ey[0] * ez[0]);
    //210
    T[5] = (ex[0] * ey[1] * ez[2] +     // 6
            ex[0] * ey[2] * ez[1] +
            ex[1] * ey[0] * ez[2] +
            ex[2] * ey[0] * ez[1] +
            ex[1] * ey[2] * ez[0] +
            ex[2] * ey[1] * ez[0]
           );
    //211
    T[6] = 2.0 * (ex[1] * ey[1] * ez[2] + // 3
                  ex[1] * ey[2] * ez[1] +
                  ex[2] * ey[1] * ez[1]
                 );
    //220
    T[7] = 2.0 * (ex[2] * ey[2] * ez[0] + // 3
                  ex[0] * ey[2] * ez[2] +
                  ex[2] * ey[0] * ez[2]);
    //221
    T[8] = 2.0 * (ex[2] * ey[2] * ez[1] + // 3
                  ex[1] * ey[2] * ez[2] +
                  ex[2] * ey[1] * ez[2]);
    //222
    T[9] = 6.0 * (ex[2] * ey[2] * ez[2]); // 1

    T *= (parity / sqrt6);

    //T/=sqrt6;
}

Particle::Particle():
    x(4), ex(3), ey(3), ez(3), Qx(6), Qy(6), Qz(6), T(10), R(3), neighbors(0)
{
    //,e2p(coord_num){
    parity = 1;
    energy = 0;
}

#include "PRE79StandardProperties.h"


void PRE79StandardProperties::CalculateMeanEnergy()
{
    double E = 0.0;
    for(int site = 0; site < lat->GetN(); site++)
    {
        E += lat->GetParticles()[site].GetEnergy();
    }
    E /= double(lat->GetN());
    energy[acc_idx] = E;
}


void PRE79StandardProperties::CalculateMeanParity()
{
    double p = 0.0;
    for(int site = 0; site < lat->GetN(); site++)
        p += lat->GetParticles()[site].GetParity();
    p /= double(lat->GetN());
    parity[acc_idx] = p;
}

void PRE79StandardProperties::CalculateMeanTensors()
{
    vect mqx(0.0, 6);
    vect mqy(0.0, 6);
    vect mqz(0.0, 6);
    vect mt(0.0, 10);
    for(int i = 0; i < lat->GetN(); i++)
    {
        mqx += lat->GetParticles()[i].GetQX();
        mqy += lat->GetParticles()[i].GetQY();
        mqz += lat->GetParticles()[i].GetQZ();
        mt += lat->GetParticles()[i].GetT();
    }
    MeanQxTensor[acc_idx] = mqx / double(lat->GetN());
    MeanQyTensor[acc_idx] = mqy / double(lat->GetN());
    MeanQzTensor[acc_idx] = mqz / double(lat->GetN());

    vect t20z(0.0, 6);
    vect t22z(0.0, 6);

    mqz /= double(lat->GetN());
    mqx /= double(lat->GetN());
    mqy /= double(lat->GetN());
    mt /= double(lat->GetN());

    t20z = std::sqrt(3. / 2.) * (mqz - Identity(3));
    t22z = std::sqrt(1. / 2.) * (mqx - mqy);
    T20T20z[acc_idx] = MatrixDotProduct(t20z, t20z);
    T22T22z[acc_idx] = MatrixDotProduct(t22z, t22z);

    vect t20x(0.0, 6);
    vect t22x(0.0, 6);

    t20x = std::sqrt(3. / 2.) * (mqx - Identity(3));
    t22x = std::sqrt(1. / 2.) * (mqy - mqz);
    T20T20x[acc_idx] = MatrixDotProduct(t20x, t20x);
    T22T22x[acc_idx] = MatrixDotProduct(t22x, t22x);

    vect t20y(0.0, 6);
    vect t22y(0.0, 6);

    t20y = std::sqrt(3. / 2.) * (mqy - Identity(3));
    t22y = std::sqrt(1. / 2.) * (mqz - mqx);
    T20T20y[acc_idx] = MatrixDotProduct(t20y, t20y);
    T22T22y[acc_idx] = MatrixDotProduct(t22y, t22y);

    T32T32[acc_idx] = Rank3Contraction(mt, mt);
}

PRE79StandardProperties::PRE79StandardProperties(): readonly(true), index(0), acc_idx(0) {}

PRE79StandardProperties::PRE79StandardProperties(const Lattice *l, int _ncycles):
    lat(l), readonly(false), ncycles(_ncycles),
    d200corz(l, _ncycles),
    d222corz(l, _ncycles),
    d220corz(l, _ncycles),
    d200corx(l, _ncycles),
    d222corx(l, _ncycles),
    d220corx(l, _ncycles),
    d200cory(l, _ncycles),
    d222cory(l, _ncycles),
    d220cory(l, _ncycles),
    d322cor(l, _ncycles),
    MeanQxTensor(_ncycles),
    MeanQyTensor(_ncycles),
    MeanQzTensor(_ncycles),
    paritycor(l, _ncycles)
{
    acc_idx = -1;
    index = 0;
    cur_b_val = 0.0;
    cur_b_idx = 0.0;
    energy.resize(ncycles, 0.0);
    parity.resize(ncycles, 0.0);
    T32T32.resize(ncycles, 0.0);

    T20T20z.resize(ncycles, 0.0);
    T22T22z.resize(ncycles, 0.0);

    T20T20x.resize(ncycles, 0.0);
    T22T22x.resize(ncycles, 0.0);

    T20T20y.resize(ncycles, 0.0);
    T22T22y.resize(ncycles, 0.0);

    b_idx.resize(ncycles, 0.0);
    b_val.resize(ncycles, 0.0);

    for(int i = 0; i < ncycles; i++)
    {
        MeanQxTensor[i].resize(6, 0.0);
        MeanQyTensor[i].resize(6, 0.0);
        MeanQzTensor[i].resize(6, 0.0);
    }

    autocorrelation_time.resize(ncycles, 0.0);
    //specific_heat.resize(ncycles);
}

void PRE79StandardProperties::Update(int idx, const StandardL2Hamiltonian *H, const AutoCorrelationTimeCalculator *ac)
{
    if(readonly) return;
    if((acc_idx + 1) >= ncycles) return;
    index = idx;

    acc_idx++;

    b_val[acc_idx] = cur_b_val;
    b_idx[acc_idx] = cur_b_idx;

    if(ac != NULL)
        autocorrelation_time[acc_idx] = ac->GetT();

    d200corz.Update();
    d220corz.Update();
    d222corz.Update();
    d200corx.Update();
    d220corx.Update();
    d222corx.Update();
    d200cory.Update();
    d220cory.Update();
    d222cory.Update();

    d322cor.Update();
    paritycor.Update();

    CalculateMeanEnergy();
    CalculateMeanParity();
    CalculateMeanTensors();

    if((acc_idx + 1) >= ncycles)
        CalculateSpecificHeat();
}

void PRE79StandardProperties::Append(const PRE79StandardProperties *_p)
{
    const PRE79StandardProperties & p = *_p;
    ncycles += p.ncycles;
    acc_idx += p.ncycles;
    index += p.ncycles;
    d200corz.Append(p.d200corz);
    d220corz.Append(p.d220corz);
    d222corz.Append(p.d222corz);
    d200corx.Append(p.d200corx);
    d220corx.Append(p.d220corx);
    d222corx.Append(p.d222corx);
    d200cory.Append(p.d200cory);
    d220cory.Append(p.d220cory);
    d222cory.Append(p.d222cory);
    d322cor.Append(p.d322cor);
    paritycor.Append(p.paritycor);

    MeanQxTensor.insert(MeanQxTensor.end(), p.MeanQxTensor.begin(), p.MeanQxTensor.end());
    MeanQyTensor.insert(MeanQyTensor.end(), p.MeanQyTensor.begin(), p.MeanQyTensor.end());
    MeanQzTensor.insert(MeanQzTensor.end(), p.MeanQzTensor.begin(), p.MeanQzTensor.end());

    int oldsize = energy.size();
    vect tmpenergy(0.0, ncycles);
    vect tmpparity(0.0, ncycles);
    vect tmpt20t20z(0.0, ncycles);
    vect tmpt22t22z(0.0, ncycles);

    vect tmpt20t20x(0.0, ncycles);
    vect tmpt22t22x(0.0, ncycles);

    vect tmpt20t20y(0.0, ncycles);
    vect tmpt22t22y(0.0, ncycles);
    vect tmpt32t32(0.0, ncycles);

    vect tmpb_val(0.0, ncycles);
    vect tmpb_idx(0.0, ncycles);

    vect tmpactime(0.0, ncycles);
    /*
        vect tmpD200(0.0,ncycles);
        vect tmpD220(0.0,ncycles);
        vect tmpD202(0.0,ncycles);
        vect tmpD222(0.0,ncycles);
         */
    for(int i = 0; i < oldsize; i++)
    {
        tmpenergy[i] = energy[i];
        tmpparity[i] = parity[i];
        tmpt20t20z[i] = T20T20z[i];
        tmpt22t22z[i] = T22T22z[i];

        tmpt20t20x[i] = T20T20x[i];
        tmpt22t22x[i] = T22T22x[i];

        tmpt20t20y[i] = T20T20y[i];
        tmpt22t22y[i] = T22T22y[i];
        tmpt32t32[i] = T32T32[i];

        tmpb_val[i] = b_val[i];
        tmpb_idx[i] = b_idx[i];

        tmpactime[i] = autocorrelation_time[i];
        /*
            tmpD200[i]=D200[i];
            tmpD220[i]=D220[i];
            tmpD202[i]=D202[i];
            tmpD222[i]=D222[i];
             */
    }

    for(int i = oldsize; i < ncycles; i++)
    {
        tmpenergy[i] = p.energy[i - oldsize];
        tmpparity[i] = p.parity[i - oldsize];
        tmpt20t20x[i] = p.T20T20x[i - oldsize];
        tmpt22t22x[i] = p.T22T22x[i - oldsize];
        tmpt20t20y[i] = p.T20T20y[i - oldsize];
        tmpt22t22y[i] = p.T22T22y[i - oldsize];
        tmpt20t20z[i] = p.T20T20z[i - oldsize];
        tmpt22t22z[i] = p.T22T22z[i - oldsize];
        tmpt32t32[i] = p.T32T32[i - oldsize];
        tmpb_val[i] = p.b_val[i - oldsize];
        tmpb_idx[i] = p.b_idx[i - oldsize];
        tmpactime[i] = p.autocorrelation_time[i - oldsize];
        /*
            tmpD200[i]=p.D200[i-oldsize];
            tmpD220[i]=p.D220[i-oldsize];
            tmpD202[i]=p.D202[i-oldsize];
            tmpD222[i]=p.D222[i-oldsize];
             */
    }
    energy.resize(ncycles, 0.0);
    parity.resize(ncycles, 0.0);
    T20T20z.resize(ncycles, 0.0);
    T22T22z.resize(ncycles, 0.0);
    T20T20x.resize(ncycles, 0.0);
    T22T22x.resize(ncycles, 0.0);
    T20T20y.resize(ncycles, 0.0);
    T22T22y.resize(ncycles, 0.0);
    T32T32.resize(ncycles, 0.0);
    b_val.resize(ncycles, 0.0);
    b_idx.resize(ncycles, 0.0);
    autocorrelation_time.resize(ncycles, 0.0);
    //D200.resize(ncycles,0.0);
    //D220.resize(ncycles,0.0);
    //D202.resize(ncycles,0.0);
    //D222.resize(ncycles,0.0);
    energy = tmpenergy;
    parity = tmpparity;
    T20T20z = tmpt20t20z;
    T22T22z = tmpt22t22z;
    T20T20x = tmpt20t20x;
    T22T22x = tmpt22t22x;
    T20T20y = tmpt20t20y;
    T22T22y = tmpt22t22y;
    T32T32 = tmpt32t32;
    b_val = tmpb_val;
    b_idx = tmpb_idx;
    autocorrelation_time = tmpactime;
    //D200=tmpD200;
    //D220=tmpD220;
    //D202=tmpD202;
    //D222=tmpD222;
}

std::valarray<bool> PRE79StandardProperties::GetReplicaMask(const double &replica) const
{
    std::cout << "Indices for replica " << replica << ": " << b_idx << std::endl;
    std::valarray<bool> mask(false, b_idx.size());

    for(int i = 0; i < b_idx.size(); i++)
        mask[i] = (b_idx[i] == replica);
    return mask;
}

PRE79StandardProperties *PRE79StandardProperties::FromReplicaMask(const PRE79StandardProperties &p, const double &replica)
{
    std::valarray<bool> mask = p.GetReplicaMask(replica);
    int n = 0;
    for(int i = 0; i < mask.size(); i++) if(mask[i]) n++;
    std::cout << "Retrieved " << n << " entries out of " << p.GetNCycles() << " for replica " << replica << "\n";

    PRE79StandardProperties * ret = new PRE79StandardProperties(p.GetLattice(), n);

    ret->energy = p.energy[mask];
    ret->parity = p.parity[mask];
    ret->T20T20z = p.T20T20z[mask];
    ret->T20T20x = p.T20T20x[mask];
    ret->T20T20y = p.T20T20y[mask];
    ret->T22T22z = p.T22T22z[mask];
    ret->T22T22x = p.T22T22x[mask];
    ret->T22T22y = p.T22T22y[mask];
    ret->T32T32 = p.T32T32[mask];
    ret->b_val = p.b_val[mask];
    ret->autocorrelation_time = p.autocorrelation_time[mask];
    ret->MeanQxTensor = p.MeanQxTensor | mask;
    ret->MeanQyTensor = p.MeanQyTensor | mask;
    ret->MeanQzTensor = p.MeanQzTensor | mask;
    (SpatialCorrelationEvolution)ret->d200corz = p.d200corz[mask];
    (SpatialCorrelationEvolution)ret->d200corx = p.d200corx[mask];
    (SpatialCorrelationEvolution)ret->d200cory = p.d200cory[mask];
    (SpatialCorrelationEvolution)ret->d220corz = p.d220corz[mask];
    (SpatialCorrelationEvolution)ret->d220corx = p.d220corx[mask];
    (SpatialCorrelationEvolution)ret->d220cory = p.d220cory[mask];
    (SpatialCorrelationEvolution)ret->d222corz = p.d222corz[mask];
    (SpatialCorrelationEvolution)ret->d222corx = p.d222corx[mask];
    (SpatialCorrelationEvolution)ret->d222cory = p.d222cory[mask];
    (SpatialCorrelationEvolution)ret->d322cor = p.d322cor[mask];
    (SpatialCorrelationEvolution)ret->paritycor = p.paritycor[mask];

    ret->acc_idx = n;

    //ret->CalculateSpecificHeat();

    return ret;
}

void PRE79StandardProperties::CalculateSpecificHeat()
{
    fluctuation = CalculateFluctuation(energy, acc_idx) * Value(lat->GetN());
}


PRE79MeanProperties::PRE79MeanProperties()
{
    int size = 0;
    mean_d200corz.resize(size, 0.0);
    mean_d222corz.resize(size, 0.0);
    mean_d220corz.resize(size, 0.0);
    mean_d200corx.resize(size, 0.0);
    mean_d222corx.resize(size, 0.0);
    mean_d220corx.resize(size, 0.0);
    mean_d200cory.resize(size, 0.0);
    mean_d222cory.resize(size, 0.0);
    mean_d220cory.resize(size, 0.0);
    mean_qx.resize(6, 0.0);
    mean_qz.resize(6, 0.0);
    mean_qy.resize(6, 0.0);

    mean_d322cor.resize(size, 0.0);
    mean_paritycor.resize(size, 0.0);
}


/// DEFUNCT - calculation of irreducible tensors, fails due to poor diagonalization
/// algorightms. The only one I know which somewhat works is in Mathematica.
void PRE79MeanProperties::CalculateMeanTensors()
{

    vect mqx(0.0, 6);
    vect mqy(0.0, 6);
    vect mqz(0.0, 6);
    mqz = MeanQxTensor();
    mqx = MeanQyTensor();
    mqy = MeanQzTensor();

    //tensory w bazie cząsteczki
    vect t20(0.0, 6);
    vect t22(0.0, 6);
    t20 = std::sqrt(1.5) * (mqz - Identity(3) / 3.);
    t22 = (mqx - mqy) / std::sqrt(2.);

    double t20evectm[3][3];
    double t20evalm[3];
    double t22evectm[3][3];
    double t22evalm[3];

    double t20mat[3][3];
    double t22mat[3][3];

    t20mat[0][0] = t20[0];
    t20mat[0][1] = t20mat[1][0] = t20[1];
    t20mat[0][2] = t20mat[2][0] = t20[2];
    t20mat[1][1] = t20[3];
    t20mat[1][2] = t20mat[2][1] = t20[4];
    t20mat[2][2] = t20[5];

    t22mat[0][0] = t22[0];
    t22mat[0][1] = t22mat[1][0] = t22[1];
    t22mat[0][2] = t22mat[2][0] = t22[2];
    t22mat[1][1] = t22[3];
    t22mat[1][2] = t22mat[2][1] = t22[4];
    t22mat[2][2] = t22[5];

    //diagonalizacja
    //dsyevj3(t20mat,t20evectm,t20evalm);
    //dsyevj3(t22mat,t22evectm,t22evalm);
    //Jacobi_Cyclic_Method((double*)t20evalm,(double*)t20evectm,(double*)t20mat,3);
    //Jacobi_Cyclic_Method((double*)t22evalm,(double*)t22evectm,(double*)t22mat,3);

    //eigen_decomposition(t20mat,t20evectm,t20evalm);
    //eigen_decomposition(t22mat,t22evectm,t22evalm);


    std::vector<evs> t20es;
    std::vector<evs> t22es;
    t20es.push_back(evs(t20evalm[0], t20evectm[0]));
    t20es.push_back(evs(t20evalm[1], t20evectm[1]));
    t20es.push_back(evs(t20evalm[2], t20evectm[2]));
    t22es.push_back(evs(t22evalm[0], t22evectm[0]));
    t22es.push_back(evs(t22evalm[1], t22evectm[1]));
    t22es.push_back(evs(t22evalm[2], t22evectm[2]));
    std::sort(t20es.begin(), t20es.end());
    std::sort(t22es.begin(), t22es.end());

    //std::cout << "ev: " << t20es[0].e << " " << t20es[1].e << " " <<  t20es[2].e << std::endl;

    vect t20a(0.0, 3);
    vect t20b(0.0, 3);
    vect t20c(0.0, 3);
    vect t22a(0.0, 3);
    vect t22b(0.0, 3);
    vect t22c(0.0, 3);

    for(int i = 0; i < 3; i++)
    {
        t20c[i] = t20es[0].v[i];
        t20a[i] = t20es[1].v[i];
        t20b[i] = t20es[2].v[i];
        t22c[i] = t22es[0].v[i];
        t22a[i] = t22es[1].v[i];
        t22b[i] = t22es[2].v[i];
    }


    //iloczyny tensorowe
    vect t20aa(0.0, 6);
    vect t20bb(0.0, 6);
    vect t20cc(0.0, 6);
    vect t22aa(0.0, 6);
    vect t22bb(0.0, 6);
    vect t22cc(0.0, 6);

    t20aa[0] = t20a[0] * t20a[0];
    t20aa[1] = t20a[1] * t20a[0];
    t20aa[2] = t20a[2] * t20a[0];
    t20aa[3] = t20a[1] * t20a[1];
    t20aa[4] = t20a[2] * t20a[1];
    t20aa[5] = t20a[2] * t20a[2];

    t20bb[0] = t20b[0] * t20b[0];
    t20bb[1] = t20b[1] * t20b[0];
    t20bb[2] = t20b[2] * t20b[0];
    t20bb[3] = t20b[1] * t20b[1];
    t20bb[4] = t20b[2] * t20b[1];
    t20bb[5] = t20b[2] * t20b[2];

    t20cc[0] = t20c[0] * t20c[0];
    t20cc[1] = t20c[1] * t20c[0];
    t20cc[2] = t20c[2] * t20c[0];
    t20cc[3] = t20c[1] * t20c[1];
    t20cc[4] = t20c[2] * t20c[1];
    t20cc[5] = t20c[2] * t20c[2];


    t22aa[0] = t22a[0] * t22a[0];
    t22aa[1] = t22a[1] * t22a[0];
    t22aa[2] = t22a[2] * t22a[0];
    t22aa[3] = t22a[1] * t22a[1];
    t22aa[4] = t22a[2] * t22a[1];
    t22aa[5] = t22a[2] * t22a[2];

    t22bb[0] = t22b[0] * t22b[0];
    t22bb[1] = t22b[1] * t22b[0];
    t22bb[2] = t22b[2] * t22b[0];
    t22bb[3] = t22b[1] * t22b[1];
    t22bb[4] = t22b[2] * t22b[1];
    t22bb[5] = t22b[2] * t22b[2];

    t22cc[0] = t22c[0] * t22c[0];
    t22cc[1] = t22c[1] * t22c[0];
    t22cc[2] = t22c[2] * t22c[0];
    t22cc[3] = t22c[1] * t22c[1];
    t22cc[4] = t22c[2] * t22c[1];
    t22cc[5] = t22c[2] * t22c[2];

    //tensory w bazie direktora
    vect t20mol(0.0, 6);
    vect t22mol(0.0, 6);
    t20mol = std::sqrt(1.5) * (t20cc - Identity(3) / 3.);
    t22mol = (t22aa - t22bb) / std::sqrt(2.);

    //std::cout << "sgn = " << sgn(t20es[0].e) << std::endl;

    mean_d200 = sgn(t20es[0].e) * MatrixDotProduct(t20mol, t20);
    mean_d220 = sgn(t20es[0].e) * MatrixDotProduct(t22mol, t20);
    mean_d202 = sgn(t20es[0].e) * MatrixDotProduct(t20mol, t22);
    mean_d222 = sgn(t20es[0].e) * MatrixDotProduct(t22mol, t22);

}

PRE79MeanProperties::PRE79MeanProperties(PRE79StandardProperties *_prop, StandardL2Hamiltonian *_H)
{
    PRE79StandardProperties & prop = *_prop;
    StandardL2Hamiltonian & H = *_H;
    mean_d200corz.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d222corz.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d220corz.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d200corx.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d222corx.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d220corx.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d200cory.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d222cory.resize(prop.GetMaxCorrLen(), 0.0);
    mean_d220cory.resize(prop.GetMaxCorrLen(), 0.0);
    mean_qx.resize(6, 0.0);
    mean_qz.resize(6, 0.0);
    mean_qy.resize(6, 0.0);

    mean_d322cor.resize(prop.GetMaxCorrLen(), 0.0);
    mean_paritycor.resize(prop.GetMaxCorrLen(), 0.0);


    energy = prop.TemporalMeanEnergyPerMolecule();
    parity = prop.TemporalMeanParity();
    parity_sus = prop.ParitySusceptibility();
    //specific_heat = prop.SpecificHeat();
    fluctuation = prop.Fluctuation();

    d200z_from_correlation = prop.Delta200ZByCorrelation();
    d222z_from_correlation = prop.Delta222ZByCorrelation();
    d200x_from_correlation = prop.Delta200XByCorrelation();
    d222x_from_correlation = prop.Delta222XByCorrelation();
    d200y_from_correlation = prop.Delta200YByCorrelation();
    d222y_from_correlation = prop.Delta222YByCorrelation();

    d200z_from_correlation_sus = prop.Delta200ZByCorrelationSusceptibility();
    d222z_from_correlation_sus = prop.Delta222ZByCorrelationSusceptibility();
    d200x_from_correlation_sus = prop.Delta200XByCorrelationSusceptibility();
    d222x_from_correlation_sus = prop.Delta222XByCorrelationSusceptibility();
    d200y_from_correlation_sus = prop.Delta200YByCorrelationSusceptibility();
    d222y_from_correlation_sus = prop.Delta222YByCorrelationSusceptibility();

    d322_from_correlation = prop.Delta322ByCorrelation();
    parity_from_correlation = prop.ParityByCorrelation();

    d322_from_correlation_sus = prop.Delta322ByCorrelationSusceptibility();
    parity_from_correlation_sus = prop.ParityByCorrelationSusceptibility();

    mean_d200corz = prop.Delta200ZMeanCorrelation();
    mean_d220corz = prop.Delta220ZMeanCorrelation();
    mean_d222corz = prop.Delta222ZMeanCorrelation();
    mean_d200corx = prop.Delta200XMeanCorrelation();
    mean_d220corx = prop.Delta220XMeanCorrelation();
    mean_d222corx = prop.Delta222XMeanCorrelation();
    mean_d200cory = prop.Delta200YMeanCorrelation();
    mean_d220cory = prop.Delta220YMeanCorrelation();
    mean_d222cory = prop.Delta222YMeanCorrelation();

    T20T20z = prop.MeanT20T20Z();
    T20T20z_sus = prop.MeanT20T20ZSusceptibility();
    T22T22z = prop.MeanT22T22Z();
    T22T22z_sus = prop.MeanT22T22ZSusceptibility();
    T20T20x = prop.MeanT20T20X();
    T20T20x_sus = prop.MeanT20T20XSusceptibility();
    T22T22x = prop.MeanT22T22X();
    T22T22x_sus = prop.MeanT22T22XSusceptibility();
    T20T20y = prop.MeanT20T20Y();
    T20T20y_sus = prop.MeanT20T20YSusceptibility();
    T22T22y = prop.MeanT22T22Y();
    T22T22y_sus = prop.MeanT22T22YSusceptibility();
    T32T32 = prop.MeanT32T32();
    T32T32_sus = prop.MeanT32T32Susceptibility();

    mean_autocorrelation_time = prop.MeanAutocorrelationTime();
    //mean_d200 = prop.MeanDelta200();
    //mean_d220 = prop.MeanDelta220();
    //mean_d202 = prop.MeanDelta202();
    //mean_d222 = prop.MeanDelta222();


    mean_d322cor = prop.Delta322MeanCorrelation();
    mean_paritycor = prop.ParityMeanCorrelation();

    mean_qx = prop.GetMeanQxTensor();
    mean_qy = prop.GetMeanQyTensor();
    mean_qz = prop.GetMeanQzTensor();

    //std::cout << mean_qx << std::endl;

    temperature = H.GetTemperature();
    tau = H.GetTau();
    lambda = H.GetLambda();
    h = H.GetH();
    kappa = H.GetKappa();
}

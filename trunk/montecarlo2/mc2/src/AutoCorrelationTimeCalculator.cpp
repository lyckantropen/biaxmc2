#include "AutoCorrelationTimeCalculator.h"


Value AutoCorrelationTimeCalculator::CalculateMeanEnergy()
{
    vect e(0.0, lattice->GetN());
    for(int site = 0; site < lattice->GetN(); site++)
    {
        e[site] = lattice->GetParticles()[site].GetEnergy();
    }
    return Mean(e, 0, 0, e.size());
}

AutoCorrelationTimeCalculator::AutoCorrelationTimeCalculator(const Lattice *l, int _nc, int _nh): lattice(l), ncycles(_nc), t(_nh / 2)
{
    ec.resize(_nc, 0.0);
    R.resize(_nh, 0.0);
    acc_idx = 0;
    //SetFile("autocorrelation");
    //SetStream(&std::cout);
}

AutoCorrelationTimeCalculator::AutoCorrelationTimeCalculator(const AutoCorrelationTimeCalculator &s)
{
    ec.resize(s.ec.size(), 0.0);
    R.resize(s.R.size(), 0.0);
    ec = s.ec;
    R = s.R;
    t = s.t;
    lattice = s.lattice;
    ncycles = s.ncycles;
    acc_idx = s.acc_idx;
}

const AutoCorrelationTimeCalculator &AutoCorrelationTimeCalculator::operator=(const AutoCorrelationTimeCalculator &s)
{
    ec.resize(s.ec.size(), 0.0);
    R.resize(s.R.size(), 0.0);
    ec = s.ec;
    R = s.R;
    t = s.t;
    lattice = s.lattice;
    ncycles = s.ncycles;
    acc_idx = s.acc_idx;
    return *this;
}

void AutoCorrelationTimeCalculator::Update()
{
    ec[acc_idx] = double(CalculateMeanEnergy());

    if(acc_idx >= (ncycles - 1))
    {
        Log() << "Autocorrelation function:\n" ;

        const int N = (ncycles - R.size() + 1);
        double em = Mean(ec);
        double emv = std::pow(Mean(ec, 0, 0, 10).Error(), 2);
        for(int k = 0; k < R.size(); k++)
        {
            double q1 = 0, q2 = 0, q3 = 0;
            for(int i = 0; i < N; i++)
            {
                q1 += ec[i] / N;
                q2 += ec[i + k] / N;
                q3 += (ec[i] - em) * (ec[i + k] - em) / N;
            }
            R[k] = q3 / emv;
            R /= R[0];
            //std::cout << q1 << ", " << q2 << ", " << q3 << std::endl;
            Log() << R[k] << std::endl;
        }
        t = R.size();
        for(int i = 0; i < R.size(); i++)
        {
            if(std::abs(R[i]) < 0.5)
            {
                t = i;
                break;
            }
        }

        Log() << "Autocorrelation time: " << t << std::endl;

        acc_idx = 0;
        return;
    }
    acc_idx++;
}

const int &AutoCorrelationTimeCalculator::GetT() const
{
    return t;
}

const vect &AutoCorrelationTimeCalculator::GetR() const
{
    return R;
}

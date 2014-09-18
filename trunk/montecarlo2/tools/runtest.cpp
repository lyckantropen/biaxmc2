#include "SimulationDB.h"
#include "Lattice.h"
#include "Metropolis.h"

int main(int argc, char ** argv)
{

    ILoggable log;
    log.SetStream(&std::cout);

    Settings set;
    set.sqlite.dir = "database.db.dir";
    set.sqlite.file = "database.db";


    log.Log() << set.ListOptions();

    set.lattice.H = 16;
    set.lattice.W = 16;
    set.lattice.L = 16;
    set.initial.isotropic = true;

    Lattice lat1(set);

    set.initial.isotropic = false;
    set.initial.biaxial = true;

    Lattice lat2(set);

    StandardL2Hamiltonian H(1.0,0.3,1,0,0);
    Metropolis met(set,&H);

    SimulationDB db(set);
    log.Log() << dynamic_cast<std::stringstream&>(db.GetDB().log()).str();
    db.GetDB().SetStream(&std::cout);


    for(int i=0; i<100; i++)
    {
        log.Log() << "this is a test: " << i << std::endl;
        int a1,b1,a2,b2;
        //lat1.Sweep(&met,a1,b1);
        db.StoreLattice(set,lat1,i);
        //lat2.Sweep(&met,a2,b2);
        db.StoreLattice(set,lat2,i);
    }
    std::cout << "FINISHED\n";
    return 0;
}

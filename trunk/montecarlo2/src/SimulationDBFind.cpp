#include "SimulationDBFind.h"


Lattice FindFinalState(const Settings &settings, bool &success)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<Lattice>    whatwegot = db.get<Lattice>(boostbase::where
                                                        (sdb.type_label, sdb.final_lattice_kw)
                                                        (sdb.H_label, settings.lattice.H)
                                                        (sdb.W_label, settings.lattice.W)
                                                        (sdb.L_label, settings.lattice.L)
                                                        (sdb.temperature_label, settings.hamiltonian.temperature)
                                                        (sdb.lambda_label, settings.hamiltonian.lambda)
                                                        (sdb.tau_label, settings.hamiltonian.tau)
                                                        (sdb.h_label, settings.hamiltonian.h)
                                                        (sdb.kappa_label, settings.hamiltonian.kappa)
                                                        //            (sdb.id_kw,settings.project.name)

                                                        );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        success = false;
        return Lattice();
    }
    else
    {
        success = true;
        return whatwegot[0];
    }
}


Lattice FindState(const Settings &settings, int production_cycle, bool &success)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<Lattice>    whatwegot = db.get<Lattice>(boostbase::where
                                                        (sdb.type_label, sdb.current_lattice_kw)
                                                        (sdb.H_label, settings.lattice.H)
                                                        (sdb.W_label, settings.lattice.W)
                                                        (sdb.L_label, settings.lattice.L)
                                                        (sdb.temperature_label, settings.hamiltonian.temperature)
                                                        (sdb.lambda_label, settings.hamiltonian.lambda)
                                                        (sdb.tau_label, settings.hamiltonian.tau)
                                                        (sdb.h_label, settings.hamiltonian.h)
                                                        (sdb.kappa_label, settings.hamiltonian.kappa)
                                                        //            (sdb.id_kw,settings.project.name)
                                                        (sdb.cycle_label, production_cycle)

                                                        );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        success = false;
        return Lattice();
    }
    else
    {
        success = true;
        return whatwegot[0];
    }
}


PRE79MeanProperties FindFinalProperties(const Settings &settings, bool &success)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<PRE79MeanProperties>    whatwegot = db.get<PRE79MeanProperties>(boostbase::where
                                                                                (sdb.type_label, sdb.final_properties_kw)
                                                                                (sdb.H_label, settings.lattice.H)
                                                                                (sdb.W_label, settings.lattice.W)
                                                                                (sdb.L_label, settings.lattice.L)
                                                                                (sdb.temperature_label, settings.hamiltonian.temperature)
                                                                                (sdb.lambda_label, settings.hamiltonian.lambda)
                                                                                (sdb.tau_label, settings.hamiltonian.tau)
                                                                                (sdb.h_label, settings.hamiltonian.h)
                                                                                (sdb.kappa_label, settings.hamiltonian.kappa)
                                                                                //            (sdb.id_kw,settings.project.name)

                                                                                );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        success = false;
        return PRE79MeanProperties();
    }
    else
    {
        success = true;
        return whatwegot[0];
    }
}


PRE79MeanProperties FindThermalizationHistory(const Settings &settings, bool &success)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<PRE79MeanProperties>    whatwegot = db.get<PRE79MeanProperties>(boostbase::where
                                                                                (sdb.type_label, sdb.thermalization_history_kw)
                                                                                (sdb.H_label, settings.lattice.H)
                                                                                (sdb.W_label, settings.lattice.W)
                                                                                (sdb.L_label, settings.lattice.L)
                                                                                (sdb.temperature_label, settings.hamiltonian.temperature)
                                                                                (sdb.lambda_label, settings.hamiltonian.lambda)
                                                                                (sdb.tau_label, settings.hamiltonian.tau)
                                                                                (sdb.h_label, settings.hamiltonian.h)
                                                                                (sdb.kappa_label, settings.hamiltonian.kappa)
                                                                                //            (sdb.id_kw,settings.project.name)

                                                                                );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        success = false;
        return PRE79MeanProperties();
    }
    else
    {
        success = true;
        return whatwegot[0];
    }
}


Lattice FindLastState(const Settings &settings, bool &success, int &cycle)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<Lattice>    whatwegot = db.get<Lattice>(boostbase::where
                                                        (sdb.type_label, sdb.current_lattice_kw)
                                                        (sdb.H_label, settings.lattice.H)
                                                        (sdb.W_label, settings.lattice.W)
                                                        (sdb.L_label, settings.lattice.L)
                                                        (sdb.temperature_label, settings.hamiltonian.temperature)
                                                        (sdb.lambda_label, settings.hamiltonian.lambda)
                                                        (sdb.tau_label, settings.hamiltonian.tau)
                                                        (sdb.h_label, settings.hamiltonian.h)
                                                        (sdb.kappa_label, settings.hamiltonian.kappa)
                                                        //            (sdb.id_kw,settings.project.name)

                                                        );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        return FindFinalState(settings, success);
    }
    else
    {
        success = true;
        cycle = whatwegot.size() * settings.simulation.measure_frequency;
        return whatwegot.back();
    }
}


Lattice FindLastStateTemperatureTolerant(const Settings &settings, bool &success, int &cycle, double tolerance)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<Lattice>    whatwegot = db.get<Lattice>(boostbase::where
                                                        (sdb.type_label, sdb.current_lattice_kw)
                                                        (sdb.H_label, settings.lattice.H)
                                                        (sdb.W_label, settings.lattice.W)
                                                        (sdb.L_label, settings.lattice.L)
                                                        //            (sdb.temperature_label,settings.hamiltonian.temperature)
                                                        (sdb.lambda_label, settings.hamiltonian.lambda)
                                                        (sdb.tau_label, settings.hamiltonian.tau)
                                                        (sdb.h_label, settings.hamiltonian.h)
                                                        (sdb.kappa_label, settings.hamiltonian.kappa)
                                                        //            (sdb.id_kw,settings.project.name)
                                                        ,
                                                        boostbase::between(sdb.temperature_label,
                                                                           settings.hamiltonian.temperature - std::abs(settings.scanning.delta) * tolerance,
                                                                           settings.hamiltonian.temperature + std::abs(settings.scanning.delta) * tolerance)
                                                        );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        return FindFinalState(settings, success);
    }
    else
    {
        success = true;
        cycle = whatwegot.size() * settings.simulation.measure_frequency;
        return whatwegot.back();
    }
}


Lattice FindLastStateFieldTolerant(const Settings &settings, bool &success, int &cycle, double tolerance)
{
    SimulationDB sdb(settings);
    boostbase::base & db = sdb.GetDB();
    std::vector<Lattice>    whatwegot = db.get<Lattice>(boostbase::where
                                                        (sdb.type_label, sdb.current_lattice_kw)
                                                        (sdb.H_label, settings.lattice.H)
                                                        (sdb.W_label, settings.lattice.W)
                                                        (sdb.L_label, settings.lattice.L)
                                                        (sdb.temperature_label, settings.hamiltonian.temperature)
                                                        (sdb.lambda_label, settings.hamiltonian.lambda)
                                                        (sdb.tau_label, settings.hamiltonian.tau)
                                                        (sdb.kappa_label, settings.hamiltonian.kappa)
                                                        //            (sdb.h_label,settings.hamiltonian.h)
                                                        //            (sdb.id_kw,settings.project.name)
                                                        ,
                                                        boostbase::between(sdb.h_label,
                                                                           settings.hamiltonian.h - tolerance,
                                                                           settings.hamiltonian.h + tolerance)
                                                        );
    //std::cout << db.log().str() << std::endl;
    if(whatwegot.size() == 0)
    {
        return FindFinalState(settings, success);
    }
    else
    {
        success = true;
        cycle = whatwegot.size() * settings.simulation.measure_frequency;
        return whatwegot.back();
    }
}

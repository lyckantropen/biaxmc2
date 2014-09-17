#include "RuntimePropertiesServer.h"
#include "PRE79Scanning.h"

/**
 * @todo równoległa produkcja i równoległe symulacje
 */
RuntimePropertiesServer::RuntimePropertiesServer(const Settings &_settings, PRE79Simulation &_sim):
    sim(&_sim),settings(&_settings),end(false),
    FIFOInterface(_settings.project.name)
{
    scanning = NULL;
    if(settings->scanning.threaded_production){
        parallel_simulation = false;
        parallel_production = true;
        serial=false;
    }
    else {
        parallel_simulation = false;
        parallel_production = false;
        serial = true;
    }
}

RuntimePropertiesServer::RuntimePropertiesServer(const Settings &_settings, PRE79Scanning &_scan):
    scanning(&_scan),settings(&_settings),end(false),
    FIFOInterface(_settings.project.name)
{
    std::cout << settings->project.name << std::endl;
    sim = NULL;
    parallel_simulation = true;
    parallel_production = false;
    serial = false;
}

RuntimePropertiesServer::RuntimePropertiesServer(const RuntimePropertiesServer &s):
    settings(s.settings),
    sim(s.sim),
    scanning(s.scanning),
    end(s.end),
    serial(s.serial),
    parallel_simulation(s.parallel_simulation),
    parallel_production(s.parallel_production),
    FIFOInterface((const FIFOInterface & )(s))
{
}

void RuntimePropertiesServer::operator()() {
    while(true){
        int l = 0;
        while(l!=-1){
            std::string line;
            l = Read(line);
            boost::trim(line);
            boost::erase_all(line,"\n");
            ////ogólne
            if(line=="settings")
                Write(settings->GetInternalLog());
            if(line=="hello")
                Write("hello\n");
            if(line=="threads"){
                std::stringstream s;
                s << omp_get_thread_num() << std::endl;
                Write(s.str());
            }
            if(line=="mode"){
                if(serial) Write("serial\n");
                if(parallel_production) Write("parallel_production\n");
                if(parallel_simulation) Write("parallel_simulation\n");
            }

            //if(line.length()>0) std::cout << line << std::endl;
            ////symulacja jednowątkowa
            if(serial){
                if(line=="lattice")
                    Write(*(sim->GetLattice()));
                if(line=="properties")
                    Write(*(sim->GetProperties()));
                if(line=="thermalization_history")
                    Write(*(sim->GetThermalizationProperties()));

                if(line=="progress"){
                    std::stringstream o;
                    o << "Thermalization: " << double(sim->GetThermalization()->GetAccIdx()+1)/double(sim->GetThermalization()->GetNCycles())*100.0 << "%\n";
                    o << "Simulation: " << double(sim->GetSimulation()->GetAccIdx()+1)/double(sim->GetSimulation()->GetNCycles())*100.0 << "%\n";
                    Write(o.str());
                }
            }
            ////
            if(parallel_production){
                const int n = sim->GetNProductions();
                std::vector<std::string> s;
                boost::split(s,line,boost::is_any_of(" "));

                /*
                if(line.length()>0)
                for(std::string & c : s){
                    std::cout << "s: " << c << "\n";
                }
                */
                 
                //std::cout << std::endl;
                if(s[0]=="production_number")
                    Write(n);
                if(s[0]=="lattice" && s.size()>1)
                    Write(*(sim->GetProduction(std::atoi(s[1].c_str()))->GetLattice()));
                if(s[0]=="lattice" && s.size()==1)
                    Write(*(sim->GetLattice()));
                if(s[0]=="properties")
                    Write(*(sim->GetProduction(std::atoi(s[1].c_str()))->GetProperties()));
                if(s[0]=="thermalization_history")
                    Write(*(sim->GetThermalizationProperties()));
                if(s[0]=="progress"){
                    std::stringstream o;
                    o << n << ": Thermalization: " << double(sim->GetThermalization()->GetAccIdx()+1)/double(sim->GetThermalization()->GetNCycles())*100.0 << "%\n";

                    for(int i=0;i<n;i++)
                        o << i << ": Simulation: " << double(sim->GetProduction(i)->GetSimulation()->GetAccIdx()+1)/double(sim->GetSimulation()->GetNCycles())*100.0 << "%\n";

                    Write(o.str());
                }
            }

            if(end) return;
        }
    }
}

void RuntimePropertiesServer::Terminate() {
    end=true;
}

#include "RuntimePropertiesServer.h"
#include "PRE79Scanning.h"

/**
 * @todo równoległa produkcja i równoległe symulacje
 */
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
                Write(settings.GetLog());
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
                foreach(std::string & c,s){
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

/* 
 * File:   main.cpp
 * Author: karol
 *
 * Created on 8 grudzień 2009, 17:42
 */

#include "base.h"
#include "boost.h"
#include "std.h"
#include "SimulationDB.h"
#include "PRE79StandardProperties.h"

/*
 * 
 */
int main(int argc, char** argv)
{
    std::string                 data_type="final_properties";
    std::string dbfile="";
    std::string dbdir="";
    boostbase::tween_t_proxy    betweens;
    boostbase::pair_t_proxy     wheres;
    std::vector<std::string>    columns;
    bool sqlite_debug=false;

    for(int i=1;i<argc;i++){
        if(std::string(argv[i])=="--db"){
            dbfile=argv[i+1];
            dbdir=argv[i+2];
        }
        if(std::string(argv[i])=="--data-type")
            data_type=argv[i+1];
        if(std::string(argv[i])=="--columns"){
            std::string columns_string = argv[i+1];
            boost::split(columns,columns_string,boost::is_any_of(","));
        }
        if(std::string(argv[i])=="--between")
            betweens(std::string(argv[i+1]),std::atof(argv[i+2]),std::atof(argv[i+3]));
        if(std::string(argv[i])=="--where")
            wheres(std::string(argv[i+1]),std::string(argv[i+2]));
        if(std::string(argv[i])=="--sqlite-debug")
            sqlite_debug=true;
        if(std::string(argv[i])=="--day"){
            pt::ptime day_start = pt::time_from_string(std::string(argv[i+1])+std::string(" 00:00:00.000000"));
            pt::ptime day_end = day_start + pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--today"){
            pt::ptime day_start = pt::second_clock::local_time();
            day_start-=day_start.time_of_day();
            pt::ptime day_end = day_start + pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--days"){
            pt::ptime day_start = pt::time_from_string(std::string(argv[i+1])+std::string(" 00:00:00.000000"));
            pt::ptime day_end = pt::time_from_string(std::string(argv[i+2])+std::string(" 00:00:00.000000"));
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--yesterday"){
            pt::ptime day_end = pt::second_clock::local_time();
            day_end-=day_end.time_of_day();
            pt::ptime day_start = day_end - pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--hours"){
            pt::ptime hour_start;
            pt::ptime hour_end;

            pt::ptime day_start = pt::second_clock::local_time();
            day_start-=day_start.time_of_day();

            hour_start = day_start + pt::duration_from_string(std::string(argv[i+1]));
            hour_end = day_start + pt::duration_from_string(std::string(argv[i+2]));

            betweens(std::string("date"),pt::to_simple_string(hour_start),pt::to_simple_string(hour_end));
        }

    }
    if(dbfile=="" || dbdir==""){
        std::cout << "please specify database file and directory with --db\n";
        std::exit(1);
    }

    wheres(std::string("data_type"),data_type);
    std::cout << std::setprecision(8) ;

    boostbase::base db(dbfile,dbdir);

    if(data_type=="final_properties" || data_type=="properties") {
        std::vector<PRE79MeanProperties> whatwegot= db.get<PRE79MeanProperties>(wheres,betweens);

        std::cout << "## some of the columns may imply an additional column with error values\n";
        std::cout << "# ";
        foreach(const std::string & col,columns){
            std::cout << col << " ";
        }
        std::cout << std::endl;

        foreach(const PRE79MeanProperties & prop,whatwegot){
            foreach(const std::string & column,columns){
                if(column=="temperature")
                    std::cout << prop.Temperature() << "\t";
                if(column=="specific_heat")
                    std::cout << prop.SpecificHeat().TableForm() << "\t";
                if(column=="energy")
                    std::cout << prop.TemporalMeanEnergyPerMolecule().TableForm() << "\t";
                if(column=="uniaxial_order")
                    std::cout << prop.UniaxialOrderByCorrelation().TableForm() << "\t";
                if(column=="biaxial_order")
                    std::cout << prop.BiaxialOrderByCorrelation().TableForm() << "\t";
                if(column=="tetrahedral_order")
                    std::cout << prop.TetrahedralOrderByCorrelation().TableForm() << "\t";
                if(column=="parity_order")
                    std::cout << prop.ParityOrderByCorrelation().TableForm() << "\t";
                if(column=="tau")
                    std::cout << prop.Tau() << "\t";
                if(column=="lambda")
                    std::cout << prop.Lambda() << "\t";
                if(column=="uniaxial_correaltion")
                    std::cout << prop.UniaxialMeanCorrelation() << "\t";
                if(column=="biaxial_correaltion")
                    std::cout << prop.BiaxialMeanCorrelation() << "\t";
                if(column=="tetrahedral_correaltion")
                    std::cout << prop.TetrahedralMeanCorrelation() << "\t";
                if(column=="parity_correaltion")
                    std::cout << prop.ParityMeanCorrelation() << "\t";
                if(column=="delta220_correaltion")
                    std::cout << prop.Delta220MeanCorrelation() << "\t";
            }
            std::cout << std::endl;
        }
    }
    if(data_type=="lattice" || data_type=="final_lattice") {
        std::vector<Lattice> whatwegot= db.get<Lattice>(wheres,betweens);

        foreach(const Lattice & lat,whatwegot){
            std::cout << "## some of the values (e.g. vectors) may span over many columns\n";
            std::cout << "# ";
            foreach(const std::string & col,columns){
                std::cout << col << " ";
            }
            std::cout << std::endl;
            for(int p=0;p<lat.GetN();p++){
                foreach(const std::string & column,columns){
                    const Particle & cp = lat.GetParticles()[p];
                    if(column=="index")
                        std::cout << p << "\t";
                    if(column=="parity")
                        std::cout << cp.GetParity() << "\t";
                    if(column=="orientation")
                        std::cout << cp.GetX() << "\t";
                    if(column=="T")
                        std::cout << cp.GetT() << "\t";
                    if(column=="energy")
                        std::cout << cp.GetEnergy() << "\t";
                    if(column=="Ex")
                        std::cout << cp.GetEX() << "\t";
                    if(column=="Ey")
                        std::cout << cp.GetEY() << "\t";
                    if(column=="Ez")
                        std::cout << cp.GetEZ() << "\t";
                    if(column=="Qx")
                        std::cout << cp.GetQX() << "\t";
                    if(column=="Qy")
                        std::cout << cp.GetQY() << "\t";
                    if(column=="Qz")
                        std::cout << cp.GetQZ() << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    if(data_type=="properties_evolution"){
        //TODO: tutaj jest segfault
        std::vector<PRE79StandardProperties> whatwegot= db.get<PRE79StandardProperties>(wheres,betweens);
        
        for(int p=0;p<whatwegot.size();p++){
            const PRE79StandardProperties & prop = whatwegot[p];
            std::cout << "## some of the values (e.g. vectors) may span over many columns\n";
            std::cout << "# ";
            foreach(const std::string & col,columns){
                std::cout << col << " ";
            }
            std::cout << std::endl;
            for(int t=0;t<prop.GetNCycles();t++){
                foreach(const std::string & column,columns){
                    if(column=="time")
                        std::cout << t << "\t";
                    // BUG: energia jest źle odczytywana/zapisywana
                    if(column=="energy")
                        std::cout << prop.EnergyEvolution()[t] << "\t";
                        //std::cout << prop.TemporalMeanEnergyPerMolecule() << "\t";
                    if(column=="uniaxial_correlation")
                        std::cout << prop.UniaxialCorrelationEvolution()[t] << "\t";
                    if(column=="biaxial_correlation")
                        std::cout << prop.BiaxialCorrelationEvolution()[t] << "\t";
                    if(column=="tetrahedral_correlation")
                        std::cout << prop.TetrahedralCorrelationEvolution()[t] << "\t";
                    if(column=="delta220_correlation")
                        std::cout << prop.Delta220CorrelationEvolution()[t] << "\t";
                    if(column=="parity_correlation")
                        std::cout << prop.ParityCorrelationEvolution()[t] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    if(sqlite_debug)
        std::cout << std::endl << db.log().str() << std::endl;

    return (EXIT_SUCCESS);
}


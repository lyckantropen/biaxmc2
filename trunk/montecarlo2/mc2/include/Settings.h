/* 
 * File:   Settings.h
 * Author: karol
 *
 * Created on 1 grudzień 2009, 21:48
 */

#ifndef _SETTINGS_H
#define	_SETTINGS_H

#include "boost.h"
#include "std.h"
#include "ILoggable.h"



class Settings:protected ILoggable {
    po::variables_map           vm;
    po::options_description     desc;
    bool TextBool(const std::string & s){
            if(boost::to_lower_copy(s)==std::string("yes"))
                    return true;
            if(boost::to_lower_copy(s)==std::string("true"))
                    return true;
            if(boost::to_lower_copy(s)==std::string("on"))
                    return true;
            if(boost::to_lower_copy(s)==std::string("no"))
                    return false;
            if(boost::to_lower_copy(s)==std::string("false"))
                    return false;
            if(boost::to_lower_copy(s)==std::string("off"))
                    return false;
            return false;
    }
public:
    struct FileNotFound:public std::exception {};
public:
    /*
     * ustawienia
     */
    struct _lattice {
        int L,W,H;
        _lattice():L(0),W(0),H(0){}
    } lattice;
    struct _hamiltonian {
        double lambda, tau, temperature, h;
        _hamiltonian():
        lambda(0.0),tau(0.0),temperature(1.0),h(0.0) {}
    } hamiltonian ;
    struct _simulation {
        bool find_thermalized;  std::string v_find_thermalized;
        bool pick_up_aborted;   std::string v_pick_up_aborted;
        bool adjust_radius;     std::string v_adjust_radius;
        bool adjust_radius_thermalization;  std::string v_adjust_radius_thermalization;
        bool measure_acceptance;std::string v_measure_acceptance;
        bool calculate_time;    std::string v_calculate_time;
        double find_thermalized_temperature_tolerance;
        double find_thermalized_h_tolerance;
        long production_cycles;
        int measure_frequency;
        int measure_acceptance_frequency;
        long thermalization_cycles;
        long supplementary_thermalization_cycles;
        int radius_adjustment_frequency;
        _simulation():
        production_cycles(1000),
        measure_frequency(10),
        thermalization_cycles(1000),
        supplementary_thermalization_cycles(0),
        radius_adjustment_frequency(100),
        v_adjust_radius("no"),
        v_find_thermalized("yes"),
        find_thermalized_temperature_tolerance(0.0),
        find_thermalized_h_tolerance(0.0),
        measure_acceptance_frequency(100),
        v_pick_up_aborted("no"),
        v_measure_acceptance("no"),
        v_adjust_radius_thermalization("yes"),
        v_calculate_time("yes")
        {}
    } simulation;
    struct _sqlite {
        std::string file;
        std::string dir;
        _sqlite():
        file("mc2.db"),
        dir("mc2") {}
    } sqlite;
    struct _output {
        bool save_configuration_evolution;  std::string v_save_configuration_evolution;
        bool save_final_configuration;      std::string v_save_final_configuration;
        bool save_final_properties;         std::string v_save_final_properties;
        bool save_properties_evolution;     std::string v_save_properties_evolution;
        bool save_intermediate_states;      std::string v_save_intermediate_states;
        int intermediate_states;
        _output():
        v_save_configuration_evolution("no"),
        v_save_final_configuration("yes"),
        v_save_final_properties("yes"),
        v_save_properties_evolution("yes"),
        v_save_intermediate_states("yes"),
        intermediate_states(10)
        {}
    } output ;
    struct _initial {
        bool biaxial;       std::string v_biaxial;
        bool righthanded;   std::string v_righthanded;
        bool isotropic;     std::string v_isotropic;
        _initial():
        v_biaxial("no"),
        v_righthanded("no"),
        v_isotropic("no")
        {}
    } initial;
    struct _scanning {
        bool enabled; std::string v_enabled;
        std::string variable;
        double start;
        double end;
        double delta;
        bool reuse_thermalized; std::string v_reuse_thermalized;
        bool threaded; std::string v_threaded;
        bool threaded_production; std::string v_threaded_production;
        bool continue_if_results_exist; std::string v_continue_if_results_exist;
        
        _scanning():
        v_reuse_thermalized("yes"),
        v_enabled("no"),
        v_threaded("no"),
        v_threaded_production("no"),
        v_continue_if_results_exist("yes")
        
        {}
    } scanning ;
    struct _openmp {
        bool dynamic; std::string v_dynamic;
        int number_of_threads;
        _openmp():
        number_of_threads(1),
        v_dynamic("no") {}
    } openmp ;
    struct _project {
        std::string name_format;
        std::string name;
    } project;
private:
    void SetupDescription(){
        desc.add_options()
        ("lattice.L",po::value<int>(&lattice.L),"Lattice longitude")
        ("lattice.W",po::value<int>(&lattice.W),"Lattice width")
        ("lattice.H",po::value<int>(&lattice.H),"Lattice height")
        ("hamiltonian.lambda",po::value<double>(&hamiltonian.lambda),"Lambda coupling constant")
        ("hamiltonian.tau",po::value<double>(&hamiltonian.tau),"Tau coupling constant")
        ("hamiltonian.temperature",po::value<double>(&hamiltonian.temperature),"Temperature")
        ("hamiltonian.h",po::value<double>(&hamiltonian.h),"Field coupling to second order tensor (along z)")
        ("simulation.production_cycles",po::value<long>(&simulation.production_cycles),"Total number of production cycles")
        ("simulation.measure_frequency",po::value<int>(&simulation.measure_frequency),"Number of cycles to skip between measurements. Must be non-zero.")
        ("simulation.thermalization_cycles",po::value<long>(&simulation.thermalization_cycles),"Number of thermalization cycles")
        ("simulation.supplementary_thermalization_cycles",po::value<long>(&simulation.supplementary_thermalization_cycles),"Number of thermalization cycles when reusing thermalized state")
        ("simulation.radius_adjustment_frequency",po::value<int>(&simulation.radius_adjustment_frequency),"Number of cycles to skip between radius adjustments. Must be non-zero.")
        ("simulation.adjust_radius",po::value<std::string>(&simulation.v_adjust_radius),"Whether to adjust MC radius during simulation")
        ("simulation.adjust_radius_thermalization",po::value<std::string>(&simulation.v_adjust_radius_thermalization),"Whether to adjust MC radius during thermalization")
        ("simulation.measure_acceptance",po::value<std::string>(&simulation.v_measure_acceptance),"Whether to measure MC acceptance rate during simulation")
        ("simulation.measure_acceptance_frequency",po::value<int>(&simulation.measure_acceptance_frequency),"Cycles to skip between measuring MC acceptance rate")
        ("simulation.find_thermalized",po::value<std::string>(&simulation.v_find_thermalized),"(yes/no) Find an already thermalized state in the database")
        ("simulation.find_thermalized_temperature_tolerance",po::value<double>(&simulation.find_thermalized_temperature_tolerance),"Temperature tolerance for thermalized state in units of scanning.delta")
        ("simulation.find_thermalized_h_tolerance",po::value<double>(&simulation.find_thermalized_h_tolerance),"Field tolerance for thermalized state in units of h")
        ("simulation.pick_up_aborted",po::value<std::string>(&simulation.v_pick_up_aborted),"(yest/no) Continue from last saved state from an aborted simulation. Relevant only when saving configurations. [DON'T USE YET]")
        ("simulation.calculate_time",po::value<std::string>(&simulation.v_calculate_time),"(yes/no) Calculate remaining simulation time")
        ("sqlite.file",po::value<std::string>(&sqlite.file),"Database file")
        ("sqlite.dir",po::value<std::string>(&sqlite.dir),"Database directory")
        ("output.save_configuration_evolution",po::value<std::string>(&output.v_save_configuration_evolution),"(yes/no) Save entire configuration evolution in production cycle (large db entry)")
        ("output.save_final_configuration",po::value<std::string>(&output.v_save_final_configuration),"(yes/no) Save final state of the lattice")
        ("output.save_final_properties",po::value<std::string>(&output.v_save_final_properties),"(yes/no) Save resulting (mean) properties of the system")
        ("output.save_properties_evolution",po::value<std::string>(&output.v_save_properties_evolution),"(yes/no) Save evolution of the properties of the system")
        ("output.save_intermediate_states",po::value<std::string>(&output.v_save_intermediate_states),"(yes/no) Save intermediate states for recovery")
        ("output.intermediate_states",po::value<int>(&output.intermediate_states),"Number of intermediate states to save")
        ("initial.biaxial",po::value<std::string>(&initial.v_biaxial),"(yes/no) Initial state is biaxial")
        ("initial.righthanded",po::value<std::string>(&initial.v_righthanded),"(yes/no) Initial state is righthanded")
        ("initial.isotropic",po::value<std::string>(&initial.v_isotropic),"(yes/no) Initial state is orientationally isotropic (can be chiral still)")
        ("scanning.enabled",po::value<std::string>(&scanning.v_enabled),"(yes/no) Should scanning be done")
        ("scanning.variable",po::value<std::string>(&scanning.variable),"Name of scanning variable (all options under hamiltonian apply)")
        ("scanning.start",po::value<double>(&scanning.start),"Starting value")
        ("scanning.end",po::value<double>(&scanning.end),"End value")
        ("scanning.delta",po::value<double>(&scanning.delta),"Interval")
        ("scanning.reuse_thermalized",po::value<std::string>(&scanning.v_reuse_thermalized),"(yes/no) Reuse last thermalized state to reduce time of simulation")
        ("scanning.threaded",po::value<std::string>(&scanning.v_threaded),"(yes/no) Threaded scanning")
        ("scanning.threaded_production",po::value<std::string>(&scanning.v_threaded_production),"(yes/no) Threaded production")
        ("scanning.continue_if_results_exist",po::value<std::string>(&scanning.v_continue_if_results_exist),"(yes/no) Skip already performed simulations and pick up last saved state")
        ("openmp.number_of_threads",po::value<int>(&openmp.number_of_threads),"Number of threads in threaded scanning")
        ("openmp.dynamic",po::value<std::string>(&openmp.v_dynamic),"(yes/no) Should the number of threads be assigned dynamically")
        ("project.name_format",po::value<std::string>(&project.name_format),"Formatted name")
        ;
    }
    void LoadBooleans(){
        output.save_configuration_evolution = TextBool(output.v_save_configuration_evolution);
        output.save_final_configuration = TextBool(output.v_save_final_configuration);
        output.save_final_properties = TextBool(output.v_save_final_properties);
        output.save_properties_evolution = TextBool(output.v_save_properties_evolution);
	output.save_intermediate_states = TextBool(output.v_save_intermediate_states);
        initial.biaxial = TextBool(initial.v_biaxial);
        initial.righthanded = TextBool(initial.v_righthanded);
        initial.isotropic = TextBool(initial.v_isotropic);
        scanning.enabled = TextBool(scanning.v_enabled);
        scanning.reuse_thermalized = TextBool(scanning.v_reuse_thermalized);
        scanning.threaded = TextBool(scanning.v_threaded);
        scanning.threaded_production = TextBool(scanning.v_threaded_production);
        scanning.continue_if_results_exist = TextBool(scanning.v_continue_if_results_exist);
        openmp.dynamic = TextBool(openmp.v_dynamic);
        simulation.find_thermalized = TextBool(simulation.v_find_thermalized);
        simulation.pick_up_aborted = TextBool(simulation.v_pick_up_aborted);
        simulation.adjust_radius = TextBool(simulation.v_adjust_radius);
        simulation.adjust_radius_thermalization = TextBool(simulation.v_adjust_radius_thermalization);
        simulation.measure_acceptance = TextBool(simulation.v_measure_acceptance);
        simulation.calculate_time = TextBool(simulation.v_calculate_time);
    }
public:
    Settings(const fs::path & file){
        if(!fs::exists(file))   throw FileNotFound();
        SetupDescription();
        std::ifstream f(file.string().c_str());
        try{
            po::basic_parsed_options<char> basic_options = po::parse_config_file(f,desc);
            po::store(basic_options,vm);
            po::notify(vm);
            LoadBooleans();
            Log() << "Parsed file " << file.string() << ". Options are:\n";
            foreach(po::basic_option<char> & o, basic_options.options){
                Log() << o.string_key << " = " << o.value[0] << std::endl;
            }

            //formatowanie nazwy
            project.name=project.name_format;
            for(std::vector<po::basic_option<char> >::iterator i=basic_options.options.begin();i!=basic_options.options.end();i++){
                if(vm.count(i->string_key)){
                    project.name = boost::replace_all_copy(project.name_format,std::string("&")+i->string_key,i->value[0]);
                }
            }

        }
        catch(po::unknown_option &e){
            Log() << "Unrecognized option found" << std::endl;
        }
        f.close();
    }
};



#endif	/* _SETTINGS_H */


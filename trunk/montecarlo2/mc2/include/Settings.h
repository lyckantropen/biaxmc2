/**
 * @file Settings.h
 * @brief Here the configuration file system sits
 * 
 * An interface for globally-available settings is created. It's based on the 
 * program_options library from Boost, which can read configuration files written
 * in the human-readable .ini format.
 * 
 * 
 */

#ifndef _SETTINGS_H
#define	_SETTINGS_H

#include "boost.h"
#include "std.h"
#include "valarray_external.h"
#include "ILoggable.h"


/**
 * @class Settings
 * @brief The class holding the static settings data
 * 
 * This class first reads and interprets the settings data from a file, and then
 * holds it in directly-applicable form. Generally, a class which needs some of
 * the settings will contain a copy of an instance of #Settings. There is no need
 * for more than one instance of Settings to exist in the entire program, so 
 * generally this instance will be created in #main() and then passed on as a
 * reference downward the hierarchy of subsequently created objects. This way a
 * global variable is avoided. 
 * 
 * 
 */
class Settings:public ILoggable {
    po::variables_map           vm;             ///< I no longer remember how this works, program_options is difficult to set up
    po::options_description     desc;           ///< see above
    /**
     * This function converts a string which holds a boolean value to it, 
     * such as "yes", "true" or "off" and returns a bool value.
     * 
     * @param s A positive or negative word, such as "yes", "on", "no", "off", "true", "false"
     * @return true or false. If the word is not recognized, false.
     */
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
    /**
     * An exception which is raised when the settings file could not be found.
     */
    struct FileNotFound:public std::exception {};
public:
    /**
     * Lattice settings
     */
    struct _lattice {
        int L; ///<lattice length
        int W; ///<lattice width
        int H; ///<lattice height
        ///Constructor
        _lattice():L(0),W(0),H(0){}
    } lattice;

     struct _lattice_boundary_conditions {
        bool periodic_boundary_condition_L; ///<lattice length
        std::string v_periodic_boundary_condition_L; ///<helper variable
        bool periodic_boundary_condition_W; ///<lattice width
        std::string v_periodic_boundary_condition_W; ///<helper variable
        bool periodic_boundary_condition_H; ///<lattice height
        std::string v_periodic_boundary_condition_H; ///<helper variable
        ///Constructor
        _lattice_boundary_conditions():periodic_boundary_condition_L(false),periodic_boundary_condition_W(false),periodic_boundary_condition_H(false),
                                       v_periodic_boundary_condition_L("no"),v_periodic_boundary_condition_W("no"),v_periodic_boundary_condition_H("no"){}
    } lattice_boundary_conditions ;
      /**
     * Parameters of the hamiltonian. If scanning applies, not all are in use.
     * E.g. if there is scanning of temperature, the temperature parameter is not
     * in use.
     */
    struct _hamiltonian {
        double lambda;          ///<lambda constant, biaxiality parameter
        double tau;             ///<strength of tetrahedratic coupling
        double temperature;     ///<temperature
        double h;               ///<field constant in second order coupling to the quadrupolar tensor
        double kappa;           ///<strength of the lattice cross coupling

        ///Constructor
        _hamiltonian():
        lambda(0.0),tau(0.0),temperature(1.0),h(0.0),kappa(0.0) {}
    } hamiltonian ;
    /**
     * Technical parameters for the simulation. 
     */
    struct _simulation {
        bool find_thermalized;          ///<should the simulation search the database for a previously thermalized state and reuse it? 
        std::string v_find_thermalized; ///<helper variable
        bool pick_up_aborted;           ///<not in use. There were plans for a possibility to recover from an interrupted simulation.   
        std::string v_pick_up_aborted;  ///<helper variable
        bool adjust_radius;             ///<should the radius of the random-walk be adjusted at runtime? Careful!     
        std::string v_adjust_radius;    ///<helper variable
        bool measure_acceptance;        ///<should the acceptance ratio be measured at runtime?
        std::string v_measure_acceptance;       ///<helper variable
        bool calculate_time;            ///<should the remaining time be calculated at runtime? Slows down the simulation!
        std::string v_calculate_time;   ///<helper variable
        double find_thermalized_temperature_tolerance;  ///<there is a possibility to pick up a state for a slightly different temperature. The tolerance is given in units of scanning.delta.
        double find_thermalized_h_tolerance;            ///<there is a possibility to pick up a state for a slightly different field. The tolerance is given in units of hamiltonian.h.
        double parity_flip_probability;                 ///<if we wanted to, we can change the flip-probability to something else than 0.5. But doesn't it violate detailed balance?
        double metropolis_lower_acceptance_limit;       ///<desired minimum acceptance ratio. E.g. 0.4
        double metropolis_higher_acceptance_limit;      ///<desired maximum acceptance ratio. E.g. 0.5
        double radius;          ///<the initial radius of the random-walk
        long production_cycles;         ///<total number of production cycles
        int measure_frequency;          ///<actually the inverse of frequency - number of cycles to skip after each measurement of system properties 
        int measure_acceptance_frequency;       ///<actually the inverse of frequency - number of cycles to skip after measuting the acceptance rate
        long thermalization_cycles;     ///<total number of initial thermalization cycles. If a previous state is recovered, this does nothing.
        long supplementary_thermalization_cycles;       ///<total number of supplementary thermalization cycles. Applies only when the initial state is recovered from the database.
        int radius_adjustment_frequency;        ///<actually the inverse of frequency - how many cycles should be skipped between every adjustment of the random-walk radius. Careful, can break detailed balance.
        int autocorrelation_length;             ///<size of the window used to measure autocorrelation of energy
        int autocorrelation_frequency;          ///<size of the window over which the average autocorrelation should be calculated
        ///Constructor
        _simulation():
        production_cycles(1000),
        measure_frequency(15),
        thermalization_cycles(1000),
        supplementary_thermalization_cycles(0),
        radius_adjustment_frequency(1000000),
        v_adjust_radius("no"),
        v_find_thermalized("yes"),
        find_thermalized_temperature_tolerance(0.0),
        find_thermalized_h_tolerance(0.0),
        measure_acceptance_frequency(100),
        v_pick_up_aborted("no"),
        v_measure_acceptance("no"),
        v_calculate_time("yes"),
        parity_flip_probability(0.5),
        radius(0.1),
        metropolis_lower_acceptance_limit(0.3),
        metropolis_higher_acceptance_limit(0.4),
        autocorrelation_length(20),
        autocorrelation_frequency(100)
        {}
    } simulation;
    /**
     * Database settings
     */
    struct _sqlite {
        std::string file; ///<where to place the database file
        std::string dir;        ///<a directory where to store data
        ///Constructor
        _sqlite():
        file("mc2.db"),
        dir("mc2") {}
    } sqlite;
    /**
     * Settings regarding how much results should be collected from the simulation.
     */
    struct _output {
        bool save_configuration_evolution;              ///<should the entire configuration-history, i.e. the state of the lattice at each time step be saved? This is very disk-intensive and very slow. The data volume is measured in tens of gigabytes.
        std::string v_save_configuration_evolution;     ///<helper variable
        bool save_final_configuration;                  ///<should the state of the lattice be saved after the simulation has ended?
        std::string v_save_final_configuration;         ///<helper variable
        bool save_final_properties;                     ///<should the final results (i.e. the final mean values) be saved?         
        std::string v_save_final_properties;            ///<helper variable
        bool save_properties_evolution;                 ///<should the entire history of instanteneous system properties be saved? Could be several gigabytes in volume, but useful if some values need to be calculated again or the convergence needs to be inspected.
        std::string v_save_properties_evolution;        ///<helper variable
        bool save_intermediate_states;                  ///<should some intermediate lattice states be saved?
        std::string v_save_intermediate_states;         ///<helper variable
        bool save_thermalization_properties;            ///<should the entire history of instanteneous system properties throught the thermalization be saved? Useful when inspecting convergence.
        std::string v_save_thermalization_properties;   ///<helper variable
        bool start_service;                             ///<currently unused. There were plans to start a mechanism for inquiring runtime system properties when the simulation hasn't yet finished.
        std::string v_start_service;                    ///<helper variable
        bool report_progress;                           ///<should the progress of the simulation be reported in a percentage form?
        std::string v_report_progress;                  ///<helper variable
        int intermediate_states;                        ///<how many intermediate states should be saved?
        ///Constructor
        _output():
        v_save_configuration_evolution("no"),
        v_save_final_configuration("yes"),
        v_save_final_properties("yes"),
        v_save_properties_evolution("yes"),
        v_save_intermediate_states("yes"),
        v_save_thermalization_properties("yes"),
        v_start_service("no"),
        v_report_progress("yes"),
        intermediate_states(10)
        {}
    } output ;
    /**
     * Settings regarding initial configuration of the lattice.
     * They can be combined in a sensible way, such that they don't contradict.
     */
    struct _initial {
        bool biaxial;                   ///<biaxial initial configuration (with the (C) axis parallel to (Z))
        std::string v_biaxial;          ///<helper variable
	bool biaxial_alt;               ///<biaxial initial configuration (with the (A) axis parallel to (Z))
        std::string v_biaxial_alt;      ///<helper variable
        bool righthanded;               ///<righthanded configuration (parity is +1 for every particle)
        std::string v_righthanded;      ///<helper variable
        bool isotropic;                 ///<isotropic configuration
        std::string v_isotropic;        ///<helper variable
        ///Constructor
        _initial():
        v_biaxial("no"),
	v_biaxial_alt("no"),
        v_righthanded("no"),
        v_isotropic("no")
        {}
    } initial;
    /**
     * Settings regarding scanning over a set of parameters.
     */
    struct _scanning {
        bool enabled;                   ///<should scanning be enabled
        std::string v_enabled;          ///<helper variable
        std::string variable;           ///<name of the variable to scan over
        double start;                   ///<initial value
        double end;                     ///<final value
        double delta;                   ///<step
        bool reuse_thermalized;         ///<if for a step in the scanning parameter a previously-thermalized state is found, proceed directly to production
        std::string v_reuse_thermalized;        ///<helper variable
        bool pass_on;                   ///<pass the final state from one step as an initial state to the next step
        std::string v_pass_on;          ///<helper variable
        bool threaded;                  ///<multiple simulations in parallel, the effect depends further on #threaded_production
        std::string v_threaded;         ///<helper variable
        bool threaded_production;       ///<production is split into several threads and the results are combined at the end to form the final ensamble
        std::string v_threaded_production;      ///<helper variable
        bool continue_if_results_exist;         ///<if set to true, if for the present value of the scanning parameter the final result is found in the database, this value will be skipped
        std::string v_continue_if_results_exist;        ///<helper variable
        bool parallel_tempering;        ///<enable parallel tempering (@todo: broken as of 02.2012)
        std::string v_parallel_tempering;       ///<helper variable
        int parallel_tempering_swapfq;  ///<inverse frequency - number of cycles between temperature exchange
        std::string v_values;           ///<helper variable
        vect values;                    ///<a list of values to iterate over
        bool separate_values;           ///<should the classical iteration with step #delta be replaced by iteration over the list #values?
        
        ///Constructor
        _scanning():
        v_reuse_thermalized("yes"),
        v_pass_on("yes"),
        v_enabled("no"),
        v_threaded("no"),
        v_threaded_production("no"),
        v_continue_if_results_exist("yes"),
        v_parallel_tempering("yes"),
        parallel_tempering_swapfq(10),
        v_values(""),
        separate_values(false)
        
        {}
    } scanning ;
    /**
     * Settins regarding parallel computing
     */
    struct _openmp {
        bool dynamic;   ///<should the number of threads be dynamic? Not recommended.
        std::string v_dynamic;  ///<helper variable
        int number_of_threads;  ///<number of OpenMP threads 
        ///Constructor
        _openmp():
        number_of_threads(1),
        v_dynamic("no") {}
    } openmp ;
    /**
     * Naming
     */
    struct _project {
        std::string name_format;        ///<project name, can contain placeholders which will be replaced by the value imported from the configuration file, such as &hamiltonian.temperature
        std::string name;               ///<project name converted from #name_format
    } project;
    /**
     * PBS stuff
     */
    struct _pbs {
        std::string queue;              ///<the PBS queue
        _pbs():
        queue("normal") {}
    } pbs;
private:
    
    /**
     * boost::program_options stuff
     */
    void SetupDescription(){
        desc.add_options()
        ("lattice.L",po::value<int>(&lattice.L),"Lattice longitude")
        ("lattice.W",po::value<int>(&lattice.W),"Lattice width")
        ("lattice.H",po::value<int>(&lattice.H),"Lattice height")
        ("lattice_boundary_conditions.periodic_boundary_condition_L",po::value<std::string>(&lattice_boundary_conditions.v_periodic_boundary_condition_L),"Choose periodic or free boundary condition in L i.e. x axis direction")
        ("lattice_boundary_conditions.periodic_boundary_condition_W",po::value<std::string>(&lattice_boundary_conditions.v_periodic_boundary_condition_W),"Choose periodic or free boundary condition in W i.e. y axis direction")
        ("lattice_boundary_conditions.periodic_boundary_condition_H",po::value<std::string>(&lattice_boundary_conditions.v_periodic_boundary_condition_H),"Choose periodic or free boundary condition in H i.e. z axis direction")
        ("hamiltonian.lambda",po::value<double>(&hamiltonian.lambda),"Lambda coupling constant")
        ("hamiltonian.kappa",po::value<double>(&hamiltonian.kappa),"Lattice coupling constant")
        ("hamiltonian.tau",po::value<double>(&hamiltonian.tau),"Tau coupling constant")
        ("hamiltonian.temperature",po::value<double>(&hamiltonian.temperature),"Temperature")
        ("hamiltonian.h",po::value<double>(&hamiltonian.h),"Field coupling to second order tensor (along z)")
        ("simulation.production_cycles",po::value<long>(&simulation.production_cycles),"Total number of production cycles")
        ("simulation.measure_frequency",po::value<int>(&simulation.measure_frequency),"Number of cycles to skip between measurements. Must be non-zero.")
        ("simulation.thermalization_cycles",po::value<long>(&simulation.thermalization_cycles),"Number of thermalization cycles")
        ("simulation.supplementary_thermalization_cycles",po::value<long>(&simulation.supplementary_thermalization_cycles),"Number of thermalization cycles when reusing thermalized state")
        ("simulation.radius_adjustment_frequency",po::value<int>(&simulation.radius_adjustment_frequency),"Number of cycles to skip between radius adjustments. Must be non-zero.")
        ("simulation.adjust_radius",po::value<std::string>(&simulation.v_adjust_radius),"Whether to adjust MC radius during simulation")
        ("simulation.measure_acceptance",po::value<std::string>(&simulation.v_measure_acceptance),"Whether to measure MC acceptance rate during simulation")
        ("simulation.measure_acceptance_frequency",po::value<int>(&simulation.measure_acceptance_frequency),"Cycles to skip between measuring MC acceptance rate")
        ("simulation.radius",po::value<double>(&simulation.radius),"Random walk radius (start value)")
        ("simulation.find_thermalized",po::value<std::string>(&simulation.v_find_thermalized),"(yes/no) Find an already thermalized state in the database")
        ("simulation.find_thermalized_temperature_tolerance",po::value<double>(&simulation.find_thermalized_temperature_tolerance),"Temperature tolerance for thermalized state in units of scanning.delta")
        ("simulation.find_thermalized_h_tolerance",po::value<double>(&simulation.find_thermalized_h_tolerance),"Field tolerance for thermalized state in units of h")
        ("simulation.pick_up_aborted",po::value<std::string>(&simulation.v_pick_up_aborted),"(yest/no) Continue from last saved state from an aborted simulation. Relevant only when saving configurations. [DON'T USE YET]")
        ("simulation.calculate_time",po::value<std::string>(&simulation.v_calculate_time),"(yes/no) Calculate remaining simulation time")
        ("simulation.parity_flip_probability",po::value<double>(&simulation.parity_flip_probability),"Probability of parity flip")
        ("simulation.metropolis_lower_acceptance_limit",po::value<double>(&simulation.metropolis_lower_acceptance_limit),"Metropolis lower acc level")
        ("simulation.metropolis_higher_acceptance_limit",po::value<double>(&simulation.metropolis_higher_acceptance_limit),"Metropolis higher acc level")
        ("simulation.autocorrelation_length",po::value<int>(&simulation.autocorrelation_length),"The maximum span of the autocorrelation of energy function")
        ("simulation.autocorrelation_frequency",po::value<int>(&simulation.autocorrelation_frequency),"Number of skips before every calculation of autocorrelation")
        ("sqlite.file",po::value<std::string>(&sqlite.file),"Database file")
        ("sqlite.dir",po::value<std::string>(&sqlite.dir),"Database directory")
        ("output.save_configuration_evolution",po::value<std::string>(&output.v_save_configuration_evolution),"(yes/no) Save entire configuration evolution in production cycle (large db entry)")
        ("output.save_final_configuration",po::value<std::string>(&output.v_save_final_configuration),"(yes/no) Save final state of the lattice")
        ("output.save_final_properties",po::value<std::string>(&output.v_save_final_properties),"(yes/no) Save resulting (mean) properties of the system")
        ("output.save_properties_evolution",po::value<std::string>(&output.v_save_properties_evolution),"(yes/no) Save evolution of the properties of the system")
        ("output.save_intermediate_states",po::value<std::string>(&output.v_save_intermediate_states),"(yes/no) Save intermediate states for recovery")
        ("output.start_service",po::value<std::string>(&output.v_start_service),"(yes/no) Start convenience service for remote inference at runtime")
        ("output.report_progress",po::value<std::string>(&output.v_report_progress),"(yes/no) Report progress")
        ("output.save_thermalization_properties",po::value<std::string>(&output.v_save_thermalization_properties),"(yes/no) Save thermalization final properties")
        ("output.intermediate_states",po::value<int>(&output.intermediate_states),"Number of intermediate states to save")
        ("initial.biaxial",po::value<std::string>(&initial.v_biaxial),"(yes/no) Initial state is biaxial")
        ("initial.biaxial_alt",po::value<std::string>(&initial.v_biaxial_alt),"(yes/no) Initial state is biaxial with b||z")
        ("initial.righthanded",po::value<std::string>(&initial.v_righthanded),"(yes/no) Initial state is righthanded")
        ("initial.isotropic",po::value<std::string>(&initial.v_isotropic),"(yes/no) Initial state is orientationally isotropic (can be chiral still)")
        ("scanning.enabled",po::value<std::string>(&scanning.v_enabled),"(yes/no) Should scanning be done")
        ("scanning.variable",po::value<std::string>(&scanning.variable),"Name of scanning variable (all options under hamiltonian apply)")
        ("scanning.start",po::value<double>(&scanning.start),"Starting value")
        ("scanning.end",po::value<double>(&scanning.end),"End value")
        ("scanning.delta",po::value<double>(&scanning.delta),"Interval")
        ("scanning.reuse_thermalized",po::value<std::string>(&scanning.v_reuse_thermalized),"DEPRECATED use scanning.pass_on")
        ("scanning.values",po::value<std::string>(&scanning.v_values),"Isolated temperature values")

        ("scanning.pass_on",po::value<std::string>(&scanning.v_pass_on),"(yes/no) When scanning with replica parallelization, pass the final state as the next initial state")

        ("scanning.threaded",po::value<std::string>(&scanning.v_threaded),"(yes/no) Threaded scanning")
        ("scanning.threaded_production",po::value<std::string>(&scanning.v_threaded_production),"(yes/no) Threaded production")
        ("scanning.continue_if_results_exist",po::value<std::string>(&scanning.v_continue_if_results_exist),"(yes/no) Skip already performed simulations and pick up last saved state")
        ("scanning.parallel_tempering",po::value<std::string>(&scanning.v_parallel_tempering),"(yes/no) Parallel tempering")
        ("scanning.parallel_tempering_swapfq",po::value<int>(&scanning.parallel_tempering_swapfq),"Swap frequency")
        ("openmp.number_of_threads",po::value<int>(&openmp.number_of_threads),"Number of threads in threaded scanning")
        ("openmp.dynamic",po::value<std::string>(&openmp.v_dynamic),"(yes/no) Should the number of threads be assigned dynamically")
        ("project.name_format",po::value<std::string>(&project.name_format),"Formatted name")
        ("pbs.queue",po::value<std::string>(&pbs.queue),"PBS queue name")
        ;
    }
    /**
     * Convert imported string values ("yes", "no") to booleans
     */
    void LoadBooleans(){
        lattice_boundary_conditions.periodic_boundary_condition_L=TextBool(lattice_boundary_conditions.v_periodic_boundary_condition_L);
        lattice_boundary_conditions.periodic_boundary_condition_W=TextBool(lattice_boundary_conditions.v_periodic_boundary_condition_W);
        lattice_boundary_conditions.periodic_boundary_condition_H=TextBool(lattice_boundary_conditions.v_periodic_boundary_condition_H);
        output.save_configuration_evolution = TextBool(output.v_save_configuration_evolution);
        output.save_final_configuration = TextBool(output.v_save_final_configuration);
        output.save_final_properties = TextBool(output.v_save_final_properties);
        output.save_properties_evolution = TextBool(output.v_save_properties_evolution);
	output.save_intermediate_states = TextBool(output.v_save_intermediate_states);
        output.save_thermalization_properties = TextBool(output.v_save_thermalization_properties);
        output.start_service = TextBool(output.v_start_service);
        output.report_progress = TextBool(output.v_report_progress);
        initial.biaxial = TextBool(initial.v_biaxial);
        initial.righthanded = TextBool(initial.v_righthanded);
        initial.isotropic = TextBool(initial.v_isotropic);
	initial.biaxial_alt = TextBool(initial.v_biaxial_alt);
        scanning.enabled = TextBool(scanning.v_enabled);
        scanning.pass_on = TextBool(scanning.v_pass_on);
        scanning.reuse_thermalized = TextBool(scanning.v_reuse_thermalized);
        scanning.threaded = TextBool(scanning.v_threaded);
        scanning.threaded_production = TextBool(scanning.v_threaded_production);
        scanning.continue_if_results_exist = TextBool(scanning.v_continue_if_results_exist);
        scanning.parallel_tempering = TextBool(scanning.v_parallel_tempering);
        openmp.dynamic = TextBool(openmp.v_dynamic);
        simulation.find_thermalized = TextBool(simulation.v_find_thermalized);
        simulation.pick_up_aborted = TextBool(simulation.v_pick_up_aborted);
        simulation.adjust_radius = TextBool(simulation.v_adjust_radius);
        simulation.measure_acceptance = TextBool(simulation.v_measure_acceptance);
        simulation.calculate_time = TextBool(simulation.v_calculate_time);
        scanning.separate_values = bool(scanning.v_values.size());
    }
    /**
     * Import scanning values when _scanning::separate_values is set to true and
     * overrides the classic scanning.
     */
    void LoadScanningValues(){
        std::vector<std::string> names;
        boost::split(names,scanning.v_values,boost::is_any_of(","));
        scanning.values.resize(names.size(),0.0);
        for(int i=0;i<names.size();i++)
            scanning.values[i]=std::atof(names[i].c_str());
    }
public:
    /**
     * Constructor. Imports the configuration file and inteprets its contents.
     * 
     * @param file Path to the configuration file.
     */
    Settings(const fs::path & file){
        //SetFile("settings");
        
        if(!fs::exists(file))   throw FileNotFound();
        SetupDescription();
        std::ifstream f(file.string().c_str());
        try{
            po::basic_parsed_options<char> basic_options = po::parse_config_file(f,desc);
            po::store(basic_options,vm);
            po::notify(vm);
            LoadBooleans();
            LoadScanningValues();
            Log() << "Parsed file " << file.string() << ". Options are:\n";
            foreach(po::basic_option<char> & o, basic_options.options){
                Log() << o.string_key << " = " << o.value[0] << std::endl;
            }

            //formatowanie nazwy
            project.name=project.name_format;
            for(std::vector<po::basic_option<char> >::iterator i=basic_options.options.begin();i!=basic_options.options.end();i++){
                if(vm.count(i->string_key)){
                    boost::replace_all(project.name,std::string("&")+i->string_key,i->value[0]);
                }
            }

        }
        catch(po::unknown_option &e){
            /*Log()*/ std::cout << "Unrecognized option found: " << e.get_option_name() << std::endl;
                      std::exit(1);
        }
        f.close();
    }
    ///Default empty constructor, needed for initialization of copies
    Settings(){}
    ///Copy constructor
    Settings(const Settings & s){
        lattice=s.lattice;
        initial=s.initial;
        openmp=s.openmp;
        output=s.output;
        pbs=s.pbs;
        project=s.project;
        scanning=s.scanning;
        simulation=s.simulation;
        sqlite=s.sqlite;
        hamiltonian=s.hamiltonian;
    }
    ///Assignment operator
    const Settings & operator=(const Settings & s){
        lattice=s.lattice;
        initial=s.initial;
        openmp=s.openmp;
        output=s.output;
        pbs=s.pbs;
        project=s.project;
        scanning=s.scanning;
        simulation=s.simulation;
        sqlite=s.sqlite;
        hamiltonian=s.hamiltonian;
        return *this;
    }
};



#endif	/* _SETTINGS_H */


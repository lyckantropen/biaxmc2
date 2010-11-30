/* 
 * File:   RuntimePropertiesServer.h
 * Author: karol
 *
 * Created on 9 wrzesień 2010, 13:06
 */

#ifndef _RUNTIMEPROPERTIESSERVER_H
#define	_RUNTIMEPROPERTIESSERVER_H

#include "FIFOInterface.h"
#include "PRE79Simulation.h"
class PRE79Scanning;

class RuntimePropertiesServer:public FIFOInterface {
    const Settings & settings;
    PRE79Simulation * sim;
    PRE79Scanning * scanning;
    bool end;
    bool serial;
    bool parallel_simulation;
    bool parallel_production;
public:
    RuntimePropertiesServer(const Settings & _settings, PRE79Simulation & _sim):
    sim(&_sim),settings(_settings),end(false),
    FIFOInterface(_settings.project.name)
    {
        scanning = NULL;
        if(settings.scanning.threaded_production){
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
    RuntimePropertiesServer(const Settings & _settings, PRE79Scanning & _scan):
    scanning(&_scan),settings(_settings),end(false),
    FIFOInterface(_settings.project.name)
    {
        std::cout << settings.project.name << std::endl;
        sim = NULL;
        parallel_simulation = true;
        parallel_production = false;
        serial = false;
    }

    RuntimePropertiesServer(const RuntimePropertiesServer & s):
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
    void operator()();

    void Terminate() {
        end=true;
    }
};




#endif	/* _RUNTIMEPROPERTIESSERVER_H */


/*
 * File:   RuntimePropertiesServer.h
 * Author: karol
 *
 * Created on 9 wrzesień 2010, 13:06
 */

#ifndef _RUNTIMEPROPERTIESSERVER_H
#define _RUNTIMEPROPERTIESSERVER_H

#include "FIFOInterface.h"
#include "PRE79Simulation.h"
class PRE79Scanning;

/// UNUSED and UNFINISHED interface whose purpose
/// was to retrieve serialized data (e.g. state of
/// the Lattice) at runtime from outside of the
/// application upon request

class RuntimePropertiesServer: public FIFOInterface
{
    const Settings * settings;
    PRE79Simulation * sim;
    PRE79Scanning * scanning;
    bool end;
    bool serial;
    bool parallel_simulation;
    bool parallel_production;
public:
    RuntimePropertiesServer(const Settings & _settings, PRE79Simulation & _sim);
    RuntimePropertiesServer(const Settings & _settings, PRE79Scanning & _scan);

    RuntimePropertiesServer(const RuntimePropertiesServer & s);
    void operator()();

    void Terminate();
};




#endif  /* _RUNTIMEPROPERTIESSERVER_H */


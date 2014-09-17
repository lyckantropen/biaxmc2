/*
 * File:   SimulationDBFind.h
 * Author: karol
 *
 * Created on 14 grudzień 2009, 13:24
 */

#ifndef _SIMULATIONDBFIND_H
#define _SIMULATIONDBFIND_H

#include "SimulationDB.h"

Lattice FindFinalState(const Settings & settings, bool & success);

Lattice FindState(const Settings & settings, int production_cycle, bool & success);

PRE79MeanProperties FindFinalProperties(const Settings & settings, bool & success);

PRE79MeanProperties FindThermalizationHistory(const Settings & settings, bool & success);

Lattice FindLastState(const Settings & settings, bool & success, int & cycle);

Lattice FindLastStateTemperatureTolerant(const Settings & settings, bool & success, int & cycle, double tolerance = 1.0);

Lattice FindLastStateFieldTolerant(const Settings & settings, bool & success, int & cycle, double tolerance = 0.1);

#endif  /* _SIMULATIONDBFIND_H */


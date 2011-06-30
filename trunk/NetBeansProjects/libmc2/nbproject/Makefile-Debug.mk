#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/382267761/Lattice.o \
	${OBJECTDIR}/_ext/382267761/Random.o \
	${OBJECTDIR}/_ext/652521094/mother.o \
	${OBJECTDIR}/_ext/652521094/RuntimePropertiesServer.o \
	${OBJECTDIR}/_ext/382267761/4DSphereRW.o \
	${OBJECTDIR}/_ext/382267761/Contractions.o \
	${OBJECTDIR}/_ext/382267761/Particle.o \
	${OBJECTDIR}/_ext/382267761/valarray_external.o \
	${OBJECTDIR}/_ext/652521094/FIFOInterface.o \
	${OBJECTDIR}/_ext/652521094/MarsagliaRNG.o \
	${OBJECTDIR}/_ext/382267761/Statistical.o


# C Compiler Flags
CFLAGS=-std=c99 -fopenmp

# CC Compiler Flags
CCFLAGS=-fopenmp
CXXFLAGS=-fopenmp

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a

${OBJECTDIR}/_ext/382267761/Lattice.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/Lattice.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.cpp

${OBJECTDIR}/_ext/382267761/Random.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/Random.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.cpp

${OBJECTDIR}/_ext/652521094/mother.o: ../../montecarlo2/mc2/src/mother.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/652521094
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/652521094/mother.o ../../montecarlo2/mc2/src/mother.cpp

${OBJECTDIR}/_ext/652521094/RuntimePropertiesServer.o: ../../montecarlo2/mc2/src/RuntimePropertiesServer.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/652521094
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/652521094/RuntimePropertiesServer.o ../../montecarlo2/mc2/src/RuntimePropertiesServer.cpp

${OBJECTDIR}/_ext/382267761/4DSphereRW.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/4DSphereRW.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.cpp

${OBJECTDIR}/_ext/382267761/Contractions.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/Contractions.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.cpp

${OBJECTDIR}/_ext/382267761/Particle.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/Particle.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.cpp

${OBJECTDIR}/_ext/382267761/valarray_external.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/valarray_external.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.cpp

${OBJECTDIR}/_ext/652521094/FIFOInterface.o: ../../montecarlo2/mc2/src/FIFOInterface.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/652521094
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/652521094/FIFOInterface.o ../../montecarlo2/mc2/src/FIFOInterface.cpp

${OBJECTDIR}/_ext/652521094/MarsagliaRNG.o: ../../montecarlo2/mc2/src/MarsagliaRNG.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/652521094
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/652521094/MarsagliaRNG.o ../../montecarlo2/mc2/src/MarsagliaRNG.cpp

${OBJECTDIR}/_ext/382267761/Statistical.o: /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/382267761
	${RM} $@.d
	$(COMPILE.cc) -g -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/382267761/Statistical.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.cpp

# Subprojects
.build-subprojects:
	cd ../libboostbase && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a

# Subprojects
.clean-subprojects:
	cd ../libboostbase && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

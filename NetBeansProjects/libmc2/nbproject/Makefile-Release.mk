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
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/RuntimePropertiesServer.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/eig3.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.o \
	${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/FIFOInterface.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/liblibmc2.a

dist/Release/GNU-Linux-x86/liblibmc2.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${RM} dist/Release/GNU-Linux-x86/liblibmc2.a
	${AR} rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibmc2.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-Linux-x86/liblibmc2.a

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Lattice.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Statistical.cpp

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/RuntimePropertiesServer.o: nbproject/Makefile-${CND_CONF}.mk ../../montecarlo2/mc2/src/RuntimePropertiesServer.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/RuntimePropertiesServer.o ../../montecarlo2/mc2/src/RuntimePropertiesServer.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Particle.cpp

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/eig3.o: nbproject/Makefile-${CND_CONF}.mk ../../montecarlo2/mc2/src/eig3.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/eig3.o ../../montecarlo2/mc2/src/eig3.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Contractions.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/valarray_external.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/4DSphereRW.cpp

${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/FIFOInterface.o: nbproject/Makefile-${CND_CONF}.mk ../../montecarlo2/mc2/src/FIFOInterface.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/_DOTDOT/_DOTDOT/montecarlo2/mc2/src/FIFOInterface.o ../../montecarlo2/mc2/src/FIFOInterface.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/mc2/src/Random.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/liblibmc2.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

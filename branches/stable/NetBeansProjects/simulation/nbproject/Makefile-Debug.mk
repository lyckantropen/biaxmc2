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
	${OBJECTDIR}/_ext/551454436/simulation.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-fopenmp -O3
CXXFLAGS=-fopenmp -O3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../libboostbase/dist/Debug/GNU-Linux-x86/liblibboostbase.a ../libmc2/dist/Debug/GNU-Linux-x86/liblibmc2.a -lboost_date_time-mt -lboost_filesystem-mt -lboost_iostreams-mt -lboost_program_options-mt -lboost_system-mt -ldl -lboost_thread-mt

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simulation

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simulation: ../libboostbase/dist/Debug/GNU-Linux-x86/liblibboostbase.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simulation: ../libmc2/dist/Debug/GNU-Linux-x86/liblibmc2.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simulation: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simulation ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/551454436/simulation.o: ../../montecarlo2/tools/simulation.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/551454436
	${RM} $@.d
	$(COMPILE.cc) -g -DSQLITE_TEMP_STORE=3 -I../../montecarlo2/mc2/include -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/551454436/simulation.o ../../montecarlo2/tools/simulation.cpp

# Subprojects
.build-subprojects:
	cd ../libboostbase && ${MAKE}  -f Makefile CONF=Debug
	cd ../libmc2 && ${MAKE}  -f Makefile CONF=Debug
	cd ../libmc2 && ${MAKE}  -f Makefile CONF=Debug
	cd ../libboostbase && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simulation

# Subprojects
.clean-subprojects:
	cd ../libboostbase && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../libmc2 && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../libmc2 && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../libboostbase && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

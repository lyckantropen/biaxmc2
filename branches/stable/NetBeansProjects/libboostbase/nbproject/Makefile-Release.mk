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
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1481557586/sqlite3.o \
	${OBJECTDIR}/_ext/1481557586/hl_md5.o \
	${OBJECTDIR}/_ext/1481557586/hl_md5wrapper.o \
	${OBJECTDIR}/_ext/1481557586/serialhash.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libboostbase

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libboostbase: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libboostbase ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/1481557586/sqlite3.o: ../../montecarlo2/boostbase/src/sqlite3.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.c) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/sqlite3.o ../../montecarlo2/boostbase/src/sqlite3.c

${OBJECTDIR}/_ext/1481557586/hl_md5.o: ../../montecarlo2/boostbase/src/hl_md5.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/hl_md5.o ../../montecarlo2/boostbase/src/hl_md5.cpp

${OBJECTDIR}/_ext/1481557586/hl_md5wrapper.o: ../../montecarlo2/boostbase/src/hl_md5wrapper.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/hl_md5wrapper.o ../../montecarlo2/boostbase/src/hl_md5wrapper.cpp

${OBJECTDIR}/_ext/1481557586/serialhash.o: ../../montecarlo2/boostbase/src/serialhash.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/serialhash.o ../../montecarlo2/boostbase/src/serialhash.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libboostbase

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

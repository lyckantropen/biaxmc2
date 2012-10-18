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
	${OBJECTDIR}/_ext/1481557586/sqlite3.o \
	${OBJECTDIR}/_ext/1481557586/hl_md5.o \
	${OBJECTDIR}/_ext/1481557586/hl_md5wrapper.o \
	${OBJECTDIR}/_ext/1481557586/serialhash.o


# C Compiler Flags
CFLAGS=-O3

# CC Compiler Flags
CCFLAGS=-fopenmp -O3
CXXFLAGS=-fopenmp -O3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibboostbase.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibboostbase.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibboostbase.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibboostbase.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibboostbase.a

${OBJECTDIR}/_ext/1481557586/sqlite3.o: ../../montecarlo2/boostbase/src/sqlite3.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.c) -g -DSQLITE_TEMP_STORE=3 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/sqlite3.o ../../montecarlo2/boostbase/src/sqlite3.c

${OBJECTDIR}/_ext/1481557586/hl_md5.o: ../../montecarlo2/boostbase/src/hl_md5.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.cc) -g -DSQLITE_TEMP_STORE=3 -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/hl_md5.o ../../montecarlo2/boostbase/src/hl_md5.cpp

${OBJECTDIR}/_ext/1481557586/hl_md5wrapper.o: ../../montecarlo2/boostbase/src/hl_md5wrapper.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.cc) -g -DSQLITE_TEMP_STORE=3 -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/hl_md5wrapper.o ../../montecarlo2/boostbase/src/hl_md5wrapper.cpp

${OBJECTDIR}/_ext/1481557586/serialhash.o: ../../montecarlo2/boostbase/src/serialhash.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1481557586
	${RM} $@.d
	$(COMPILE.cc) -g -DSQLITE_TEMP_STORE=3 -I../../montecarlo2/boostbase/include -I../../montecarlo2/boostbase/include/hashlib -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1481557586/serialhash.o ../../montecarlo2/boostbase/src/serialhash.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibboostbase.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

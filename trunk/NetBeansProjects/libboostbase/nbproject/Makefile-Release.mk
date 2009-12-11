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
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5wrapper.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/serialhash.o \
	${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5.o

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
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/libboostbase

dist/Release/GNU-Linux-x86/libboostbase: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libboostbase ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5wrapper.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5wrapper.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5wrapper.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5wrapper.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/serialhash.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/serialhash.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/serialhash.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/serialhash.cpp

${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5.o: nbproject/Makefile-${CND_CONF}.mk /home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5.o /home/karol/NetBeansProjects/biaxmc2/montecarlo2/boostbase/src/hl_md5.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/libboostbase

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

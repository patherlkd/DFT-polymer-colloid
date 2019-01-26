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
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/colDFT.o \
	${OBJECTDIR}/colDFTfreeenergy.o \
	${OBJECTDIR}/colDFTfuncs.o \
	${OBJECTDIR}/colDFTgetters.o \
	${OBJECTDIR}/colDFTmain.o \
	${OBJECTDIR}/colDFTpolymerdensity.o \
	${OBJECTDIR}/colDFTpotentials.o \
	${OBJECTDIR}/colDFTsetters.o \
	${OBJECTDIR}/colDFTweights.o \
	${OBJECTDIR}/random.o \
	${OBJECTDIR}/useful.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-I /homes/ldavis/Eigen/
CXXFLAGS=-I /homes/ldavis/Eigen/

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dftpolymercolloid

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dftpolymercolloid: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dftpolymercolloid ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/colDFT.o: colDFT.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFT.o colDFT.cpp

${OBJECTDIR}/colDFTfreeenergy.o: colDFTfreeenergy.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTfreeenergy.o colDFTfreeenergy.cpp

${OBJECTDIR}/colDFTfuncs.o: colDFTfuncs.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTfuncs.o colDFTfuncs.cpp

${OBJECTDIR}/colDFTgetters.o: colDFTgetters.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTgetters.o colDFTgetters.cpp

${OBJECTDIR}/colDFTmain.o: colDFTmain.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTmain.o colDFTmain.cpp

${OBJECTDIR}/colDFTpolymerdensity.o: colDFTpolymerdensity.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTpolymerdensity.o colDFTpolymerdensity.cpp

${OBJECTDIR}/colDFTpotentials.o: colDFTpotentials.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTpotentials.o colDFTpotentials.cpp

${OBJECTDIR}/colDFTsetters.o: colDFTsetters.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTsetters.o colDFTsetters.cpp

${OBJECTDIR}/colDFTweights.o: colDFTweights.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTweights.o colDFTweights.cpp

${OBJECTDIR}/random.o: random.C
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/random.o random.C

${OBJECTDIR}/useful.o: useful.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/useful.o useful.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

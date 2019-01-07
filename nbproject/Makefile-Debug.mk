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
	${OBJECTDIR}/_ext/b2c86978/DFTmain.o \
	${OBJECTDIR}/_ext/b2c86978/random.o \
	${OBJECTDIR}/_ext/b2c86978/useful.o \
	${OBJECTDIR}/colDFT.o \
	${OBJECTDIR}/colDFTComps.o \
	${OBJECTDIR}/colDFTSetters.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dftpolymercolloid

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dftpolymercolloid: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dftpolymercolloid ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/_ext/b2c86978/DFTmain.o: ../../Downloads/ColloidDFT/DFTmain.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/b2c86978
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/b2c86978/DFTmain.o ../../Downloads/ColloidDFT/DFTmain.cpp

${OBJECTDIR}/_ext/b2c86978/random.o: ../../Downloads/ColloidDFT/random.C
	${MKDIR} -p ${OBJECTDIR}/_ext/b2c86978
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/b2c86978/random.o ../../Downloads/ColloidDFT/random.C

${OBJECTDIR}/_ext/b2c86978/useful.o: ../../Downloads/ColloidDFT/useful.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/b2c86978
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/b2c86978/useful.o ../../Downloads/ColloidDFT/useful.cpp

${OBJECTDIR}/colDFT.o: colDFT.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFT.o colDFT.cpp

${OBJECTDIR}/colDFTComps.o: colDFTComps.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTComps.o colDFTComps.cpp

${OBJECTDIR}/colDFTSetters.o: colDFTSetters.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/colDFTSetters.o colDFTSetters.cpp

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

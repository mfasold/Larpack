CXX=g++
CXXFLAGS=-Wall -o3 -ansi -pedantic -fopenmp
# -gstabs for debugging
AFFYSDK_DIR=AffymetrixSdk
# AFFYSDK_DIR=/scr/achterl/mario/fusion/affy/sdk
# AFFYSDK_DIR=../../../affy/sdk
AFFYSDK_INCLUDES=-I${AFFYSDK_DIR} -I${AFFYSDK_DIR}/file -I${AFFYSDK_DIR}/portability
AFFYSDK_DEPENDS=\
	${AFFYSDK_DIR}/file/CELFileData.o \
	${AFFYSDK_DIR}/file/CELFileWriter.o \
	${AFFYSDK_DIR}/file/FileWriter.o \
	${AFFYSDK_DIR}/file/BPMAPFileData.o \
	${AFFYSDK_DIR}/file/CDFFileData.o \
	${AFFYSDK_DIR}/file/FileIO.o 
# \
#	${AFFYSDK_DIR}/file/GridControlData.o \


INCLUDE=-I/usr/include/boost -I/usr/include -IMath -IMicroarrayData -IOligoMethod -IUtility -IOpticalBackgroundCorrection ${AFFYSDK_INCLUDES}
#INCLUDE=-IMath -IMicroarrayData -IOligoMethod -IUtility -IOpticalBackgroundCorrection ${AFFYSDK_INCLUDES}
LIBS=-lgsl -lgslcblas -lm # -lboost_program_options
TESTLIBS=-lboost_unit_test_framework

DEPENDS=Utility/anyoption.o Utility/StringUtil.o Utility/StlUtil.o Utility/SequenceUtil.o Utility/DataCollector.o Math/MathUtil.o  MicroarrayData/Probeset.o OligoMethod/HookModel.o OligoMethod/SensitivityProfile.o OligoMethod/SubsetSensitivityProfile.o OligoMethod/DetectKinkPoint.o MicroarrayData/ProbesetComposition.o MicroarrayData/MicroarrayExperiment.o MicroarrayData/ProbeFilter.o MicroarrayData/SnpProbe.o MicroarrayData/PmMmProbe.o MicroarrayData/Probe.o OligoMethod/HookcurveAnalyzer.o OligoMethod/HookcurveStatistics.o OligoMethod/ProgramOptions.o OligoMethod/OligoController.o OpticalBackgroundCorrection/BackgroundSubtraction.o MicroarrayData/Chip.o OligoMethod/ExpressionMeasure.o ${AFFYSDK_DEPENDS}

# recompile if dependecies have changed
arrayCorrect: arrayCorrect.o ${DEPENDS}
	${CXX} ${CXXFLAGS} ${INCLUDE} -o arrayCorrect  ${DEPENDS} arrayCorrect.o ${LIBS}

snpCorrect: snpCorrect.o ${DEPENDS}
	${CXX} ${CXXFLAGS} ${INCLUDE} -o snpCorrect  ${DEPENDS} snpCorrect.o ${LIBS}

pmOnlyCorrect: pmOnlyCorrect.o ${DEPENDS}
	${CXX} ${CXXFLAGS} ${INCLUDE} -o pmOnlyCorrect  ${DEPENDS} pmOnlyCorrect.o ${LIBS}

computeHookStats: computeHookStats.o ${DEPENDS}
	${CXX} ${CXXFLAGS} ${INCLUDE} -o computeHookStats ${DEPENDS} computeHookStats.o ${LIBS}

# Make test-cases
test: test.o ${DEPENDS}
	${CXX} ${CXXFLAGS} ${INCLUDE} -o test  ${DEPENDS} test.o ${LIBS} ${TESTLIBS}

%.o : %.cpp
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $< -o $@

clean:
	rm -f arrayCorrect arrayCorrect.o ${DEPENDS}
	rm -f snpCorrect snpCorrect.o 
	rm -f pmOnlyCorrect pmOnlyCorrect.o

doc:
	doxygen Doxyconfig

tags:
	etags `find -name '*.cpp' -or -name '*.hpp'`

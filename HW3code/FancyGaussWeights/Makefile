.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================

C_SOURCES = gaussQuadWeight.cpp
MY_INCLUDES = gaussQuadWeight.h gaussElim.h


TARGET2 = gaussQuad
C_OBJS2 = gaussQuadWeight.o


ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)


CCX = g++
CXXFLAGS = -g -Wall
#-std=c99

#============================================================
all: $(TARGET1) $(TARGET2)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

#$(TARGET1) :   $(C_OBJS1)
#	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@
$(TARGET2) :   $(C_OBJS2)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@ 

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================


clean:
	rm -f $(TARGET1) $(TARGET2) $(C_OBJS1) $(C_OBJS2) *~

tar: $(ALL_SOURCES) $(DATA_FILES)
	tar -cvf FancyGaussWeights.tar  $(ALL_SOURCES)   




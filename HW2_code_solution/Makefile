.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================

TARGET1 = search
C_OBJS1 = search.o 

TARGET2 = newton
C_OBJS2 = bisection_vs_newton.o 

TARGET3 = getZeros
C_OBJS3 = getZeros.o 


C_SOURCES =  search.cpp bisection_vs_newton.cpp getZeros.cpp
MY_INCLUDES = search.h sort.h

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

CCX = g++
CXXFLAGS = -g -Wall
SHELL:=/bin/bash
#-std=c99

#============================================================
all: 	$(TARGET1) $(TARGET2)  $(TARGET3)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

$(TARGET2) :   $(C_OBJS2)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@ 

$(TARGET3) :   $(C_OBJS3)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@ 


# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================


clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3)  $(C_OBJS1) $(C_OBJS2)  $(C_OBJS3) *~

tar: $(ALL_SOURCES)
	tar cvf HW2_code.tar $(ALL_SOURCES)




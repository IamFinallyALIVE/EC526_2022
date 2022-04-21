#============================================================
TARGET	=  cluster

SOURCES = cluster.cpp

OBJS =  cluster.o


MY_INCLUDES = rng.h

COMP = g++

FLAGS =  -I/usr/X11/include

LIB =  -framework OpenGL -framework GLUT -framework Cocoa
#============================================================
all: $(SOURCES) $(TARGET)

%.o : %.cpp 
	$(COMP)  -c $(FLAGS) $<  

%.o : %.cu $(MY_INCLUDES)
	$(COMP)  -c  $(FLAGS) $<  

$(TARGET) :  $(OBJS)
	$(COMP) $(CXXFLAGS) $(LIB) $(OBJS)  -o $@

# Implicit rules: 
# $@ = target name, $< = first prerequisite name, $^ = all prerequisites 
#============================================================

ALL_SOURCES = Makefile $(SOURCES) $(MY_INCLUDES)

NvidiaDoc = 

clean:
	rm -f $(TARGET) $(OBJS) core 

tar: $(ALL_SOURCES) $(NvidiaDoc)
	tar cvf $(TARGET).tar $(ALL_SOURCES) $(NvidiaDoc)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)



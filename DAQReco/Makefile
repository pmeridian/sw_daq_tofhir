CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)
INC = $(shell pwd)
Repo = $(shell git rev-parse --show-toplevel)
SRC = $(Repo)/src

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs)
CPPFLAGS += -g -std=c++1y

TARGET = ConvertTOFPETSinglesToEvents
SRC = ConvertTOFPETSinglesToEvents.cc
OBJ = $(SRC:.cc=.o)
head= Tofhir.h


all :  $(TARGET)


$(TARGET) : $(OBJ)
	@echo $@
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(head) $(LDFLAGS) 


%.o : %.cc
	@echo $@
	$(CXX) $(CPPFLAGS) -o $@ -c $<
clean :
	rm -f *.o $(TARGET) *~

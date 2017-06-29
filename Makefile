# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler, linker, ... flags:
CC_FLAGS := $(shell root-config --cflags)
LDFLAGS := $(shell root-config --ldflags)
LIBS := $(shell root-config --glibs) -lMathMore

START := $(shell rm -rf LinkDef.*)

CPP_FILES := $(wildcard *.cpp)

HEADERS :=  $(addprefix ./,$(notdir $(CPP_FILES:.cpp=.h)))

SRC_MAIN=main.C gaus2.C

OBJ_FILES := $(addprefix ./,$(notdir $(CPP_FILES:.cpp=.o)))
OBJ_MAIN=$(SRC_MAIN:.C=.o)

# the build target executable:
TEST = test.exe
GAUS2 = gaus2.exe

# build the executable (remove LinkDef at the end to avoid conflict in the next "make" call)
all:    clean test gaus2

# generate dictionary

LinkDef.cpp: $(HEADERS)
	rootcint -f LinkDef.cpp -c $(CC_FLAGS) $(HEADERS) My_LinkDef.h

# build reco executable

test:   $(OBJ_FILES) $(OBJ_MAIN) LinkDef.o main.o
	$(CC) $(CC_FLAGS) -o $(TEST) $(OBJ_FILES) LinkDef.o main.o $(LFLAGS) $(LIBS)

gaus2:  $(OBJ_FILES) $(OBJ_MAIN) LinkDef.o gaus2.o
	$(CC) $(CC_FLAGS) -o $(GAUS2) $(OBJ_FILES) LinkDef.o gaus2.o $(LFLAGS) $(LIBS)

#compilation commands
.cpp.o:
	g++ $(CC_FLAGS) -c -o $@ $<
.C.o:
	g++ $(CC_FLAGS) -c -o $@ $<

clean:
	rm -rf LinkDef.* *.o bin/*.exe *~





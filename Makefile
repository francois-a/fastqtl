#PLEASE SPECIFY THE R path here where you built the R math library standalone 
RMATH=/home/jjzhu/src/R-3.2.0/src
CBOOST=/home/jjzhu/src/boost_1_58_0
CNPY=/home/jjzhu/src/cnpy



#compiler
CXX=g++ -std=c++11

#internal paths
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)
PATH_TABX=lib/Tabix
PATH_EIGN=lib/Eigen

#compiler flags
CXXFLAG_OPTI=-O3 -D_FAST_CORRELATION
CXXFLAG_DEBG=-g
CXXFLAG_WARN=-Wall -Wextra -Wno-sign-compare
CXXFLAG_MACX=-mmacosx-version-min=10.7 -stdlib=libc++

#linker flags
LDFLAG_OPTI=-O3
LDFLAG_DEBG=-g
LDFLAG_MACX=-mmacosx-version-min=10.7 -stdlib=libc++

#includes
# INC_BASE=-Isrc -I$(PATH_TABX) -I$(PATH_EIGN)
INC_BASE=-Isrc -I$(PATH_TABX) -I$(PATH_EIGN) -I$(CBOOST)/include -I$(CNPY)/include
INC_MATH=-I$(RMATH)/include/
INC_MACX=-I/usr/local/include/

#libraries
#LIB_BASE=-lm -lboost_iostreams -lboost_program_options -lz -lgsl -lblas
#LIB_BASE=-lm -lz -lbz2 -lboost_iostreams -lboost_program_options -lgslcblas -lgsl -lblas
LIB_BASE=-lm -lz -lboost_iostreams -lboost_program_options -lgsl -lgslcblas -lcnpy -L$(CBOOST)/lib -L$(CNPY)/lib 
LIB_MATH=$(RMATH)/nmath/standalone/libRmath.a
LIB_TABX=$(PATH_TABX)/libtabix.a
LIB_MACX=-L/usr/local/lib/

#files (binary, objects, headers & sources)
FILE_BIN=bin/fastQTL
FILE_O=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
FILE_H=$(shell find src -name *.h)
FILE_CPP=$(shell find src -name *.cpp)

#default
all: linux

#linux release
linux: CXXFLAG=$(CXXFLAG_OPTI) $(CXXFLAG_WARN)
linux: IFLAG=$(INC_BASE) $(INC_MATH)
linux: LIB=$(LIB_MATH) $(LIB_TABX) $(LIB_BASE)
linux: LDFLAG=$(LDFLAG_OPTI)  
linux: $(FILE_BIN)

#macos release
macos: CXXFLAG=$(CXXFLAG_OPTI) $(CXXFLAG_WARN) $(CXXFLAG_MACX)
macos: IFLAG=$(INC_BASE) $(INC_MACX) $(INC_MATH)
macos: LIB=$(LIB_MACX) $(LIB_MATH) $(LIB_TABX) $(LIB_BASE)
macos: LDFLAG=$(LDFLAG_OPTI) $(LDFLAG_MACX)  
macos: $(FILE_BIN)

#debug release
debug: CXXFLAG=$(CXXFLAG_DEBG) $(CXXFLAG_WARN)
debug: IFLAG=$(INC_BASE) $(INC_MATH)
debug: LIB=$(LIB_MATH) $(LIB_TABX) $(LIB_BASE)
debug: LDFLAG=$(LDFLAG_DEBG)
debug: $(FILE_BIN)

#compilation
$(LIB_TABX):
	cd $(PATH_TABX) && make && cd ../..

$(FILE_BIN): $(FILE_O) $(LIB_TABX)
	$(CXX) $(LDFLAG) $^ $(LIB) -o $@

obj/%.o: %.cpp $(FILE_H)
	$(CXX) $(CXXFLAG) -o $@ -c $< $(IFLAG)

clean: 
	rm -f obj/*.o $(FILE_BIN)
	
cleanall: clean 
	cd $(PATH_TABX) && make clean && cd ../..

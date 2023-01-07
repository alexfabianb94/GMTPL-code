SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio128/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio128/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl

BUILDPATH	= build/
DATAEXAMPLE	= ../data/GMTPL/cplex/example-data/parametros.dat 
BINARIO	= GMTPL.out

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I./include
BOOSTFLAG = -lboost_program_options 
#-static

all: $(BINARIO) clean_aux


# ------------------------------------------------------------

clean :
	/bin/rm -rf $(BUILDPATH)$(BINARIO) $(BUILDPATH)*.dat $(BUILDPATH)*.sh
	/bin/rm -rf *.o 
	
clean_aux:
	/bin/rm -rf *.o
	
# ------------------------------------------------------------
#
# The examples
#
# 

$(BINARIO): data.o solver_flow.o solver_multi.o main.o solver_mtz.o solver_mtz_m.o solver_multi_m.o lib/libboost_program_options.a
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(BUILDPATH)$(BINARIO) $? $(CCLNFLAGS) $(BOOSTFLAG)
	
	
%.o: %.cpp
	$(CCC) -c $(CCFLAGS) $< -o $@
	
main.o : main.cpp
	$(CCC) -c $(CCFLAGS) $< -o $@ $(BOOSTFLAG)

example-data: $(DATAEXAMPLE)
	cp $(DATAEXAMPLE) $(BUILDPATH)

# Local Variables:
# mode: makefile
# End:

################################################################################
# When you adapt this makefile to compile your CPLEX programs
# set CPLEXDIR to the directories where CPLEX is installed:
#
CPLEXDIR  = ../../../../encs/pkg/cplex-12.7.1/root/cplex
#
# Link options and libraries for 64 Bits CPLEX
SYSTEM    = x86-64_linux
LIBFORMAT = static_pic
#
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXBINDIR = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXINCDIR = $(CPLEXDIR)/include/ilcplex/
#
CLNDIRS  = -L$(CPLEXLIBDIR)
CLNFLAGS = -lcplex -lm -lpthread
#
################################################################################


# Compiler selection 
CC     = gcc -O2
CFLAGS = -m64 -fPIC -fno-strict-aliasing -I$(CPLEXINCDIR)

#------------------------------------------------------------
#  make all     : to compile the examples. 
#  make execute : to compile and execute the examples.
#------------------------------------------------------------

main: main.o read_data.o output_data.o TwoIndex.o ThreeIndex.o FourIndex.o\
      Heuristics.o Benders.o Benders_master.o Benders_subproblem.o \
      CLM_utils.o CLM_graph.o CLM_OCSTP.o
	$(CC) $(CFLAGS) $(CLNDIRS) $^ -o $@ $(CLNFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

main.o:                 main.c               headers.h
read_data.o:            read_data.c          headers.h
output_data.o:          output_data.c        headers.h
TwoIndex.o:             TwoIndex.c           headers.h
ThreeIndex.o:           ThreeIndex.c         headers.h
FourIndex.o:            FourIndex.c          headers.h
Heuristics.o:           Heuristics.c         headers.h
Benders.o:              Benders.c            headers.h
Benders_master.o:       Benders_master.c     headers.h
Benders_subproblem.o:   Benders_subproblem.c headers.h
CLM_utils.o:            CLM_utils.c          headers.h
CLM_graph.o:            CLM_graph.c          headers.h
CLM_OCSTP.o:            CLM_OCSTP.c          headers.h

all:
	make main
    
run: all
	./main ./configuration.txt ./instances.txt

clean:
	/bin/rm -rf *.o *~
	/bin/rm -rf main
 

# makefile 
#*****************************************************************
#system dependent defaults
DEBUG = -g
OPTIMIZE = -O
MPICC = mpic++
GCC = g++
RM = rm -f
#CFLAGS   = $(OPTIMIZE) -Wall $(DEBUG)

# Macro definition
# [100, 500, 1000, 5000, 10000, 50000]
POPSIZE = 5000
MAXGENS = 100
CFLAGS   = -Wall $(DEBUG) -DDEFINED_BY_MAKE -DMAXGENS=$(MAXGENS) -DPOPSIZE=$(POPSIZE)
#CFLAGS   = -Wall $(DEBUG)
#*****************************************************************
# source files to be compiled and linked
SGA = sga
PGA = pga_global

PGA_EXEC = $(PGA)_Report_P$(POPSIZE)
SGA_EXEC = $(SGA)_Report_P$(POPSIZE)

#*****************************************************************
#all: $(SGA_EXEC) $(PGA_EXEC)
all: $(PGA_EXEC)
	
$(PGA_EXEC): $(PGA).o
	$(MPICC) -o $@ $^
	mv $@ bin/
	
$(PGA).o: $(PGA).cpp
	$(MPICC) $(CFLAGS) -c $<
	
$(SGA_EXEC): $(SGA).o
	$(GCC) -o $@ $^
	mv $@ bin/
	
$(SGA).o: $(SGA).cpp
	$(GCC) $(CFLAGS) -c $<
#*****************************************************************
clean:
#	$(RM) bin/$(SGA_EXEC) bin/$(PGA_EXEC)
	$(RM) bin/$(PGA_EXEC)
	find . -name "*~" -exec rm {} \;
	find . -name "*.o" -exec rm {} \;

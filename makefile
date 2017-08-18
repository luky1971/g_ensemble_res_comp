# Set CC and CXX to the compilers that were used to compile Gromacs
CFLAGS += -std=c99 -O3
# CFLAGS += -std=c99 -g -DEC_DEBUG
# CFLAGS += -std=c99 -g

PARALLEL = 1

GROMACS = /usr/local/gromacs
VGRO = 5
SVM = extern/libsvm-3.20

SRC = src
BUILD = build
INSTALL = /usr/local/bin

LIBS = -ldl

ifeq ($(VGRO),5)
INCGRO = \
	-I$(GROMACS)/include/ \
	-I$(GROMACS)/include/gromacs/utility \
	-I$(GROMACS)/include/gromacs/fileio \
	-I$(GROMACS)/include/gromacs/commandline \
	-I$(GROMACS)/include/gromacs/legacyheaders \
	-I$(GROMACS)/include/gromacs/topology
LINKGRO = -L$(GROMACS)/lib
LIBGRO = -lgromacs
DEFV5 = -D GRO_V5
else
INCGRO = -I$(GROMACS)/include/gromacs
LINKGRO = -L$(GROMACS)/lib
LIBGRO = -lgmx
endif

ifneq ($(PARALLEL),0)
CFLAGS += -fopenmp
endif

.PHONY: install clean

$(BUILD)/g_ensemble_comp: $(BUILD)/g_ensemble_comp.o $(BUILD)/ensemble_comp.o
	make svm.o -C $(SVM) \
	&& $(CXX) $(CFLAGS) -o $(BUILD)/g_ensemble_comp $(BUILD)/g_ensemble_comp.o $(BUILD)/ensemble_comp.o \
	$(SVM)/svm.o $(LINKGRO) $(LIBGRO) $(LIBS)

install: $(BUILD)/g_ensemble_comp
	install $(BUILD)/g_ensemble_comp $(INSTALL)

$(BUILD)/g_ensemble_comp.o: $(SRC)/g_ensemble_comp.c $(SRC)/ensemble_comp.h
	$(CC) $(CFLAGS) -o $(BUILD)/g_ensemble_comp.o -c $(SRC)/g_ensemble_comp.c $(DEFV5) $(INCGRO) -I$(SVM)

$(BUILD)/ensemble_comp.o: $(SRC)/ensemble_comp.c $(SRC)/ensemble_comp.h
	$(CC) $(CFLAGS) -o $(BUILD)/ensemble_comp.o -c $(SRC)/ensemble_comp.c $(DEFV5) $(INCGRO) -I$(SVM)

clean:
	make clean -C $(SVM) \
	&& rm -f $(BUILD)/*.o $(BUILD)/g_ensemble_comp

# TODO: MAKE THIS CLEANER
# Set CC and CXX to the compilers that were used to compile Gromacs
CFLAGS += -std=c99 -O3
# CFLAGS += -std=c99 -g -DEC_DEBUG
# CFLAGS += -std=c99 -g

PARALLEL = 1

GROMACS = /usr/local/gromacs
VGRO = 5
SVM = extern/libsvm-3.20

GKUT = extern/gkut

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

$(BUILD)/g_ensemble_res_comp: $(BUILD)/g_ensemble_res_comp.o $(BUILD)/ensemble_res_comp.o gkut
	make svm.o -C $(SVM) \
	&& $(CXX) $(CFLAGS) -o $(BUILD)/g_ensemble_res_comp $(BUILD)/g_ensemble_res_comp.o $(BUILD)/ensemble_res_comp.o \
	$(SVM)/svm.o $(GKUT)/build/gkut_io.o $(GKUT)/build/gkut_log.o $(LINKGRO) $(LIBGRO) $(LIBS)

install: $(BUILD)/g_ensemble_res_comp
	install $(BUILD)/g_ensemble_res_comp $(INSTALL)

$(BUILD)/g_ensemble_res_comp.o: $(SRC)/g_ensemble_res_comp.c $(SRC)/ensemble_res_comp.h
	$(CC) $(CFLAGS) -o $(BUILD)/g_ensemble_res_comp.o -c $(SRC)/g_ensemble_res_comp.c $(DEFV5) $(INCGRO) -I$(SVM) -I$(GKUT)/include

$(BUILD)/ensemble_res_comp.o: $(SRC)/ensemble_res_comp.c $(SRC)/ensemble_res_comp.h
	$(CC) $(CFLAGS) -o $(BUILD)/ensemble_res_comp.o -c $(SRC)/ensemble_res_comp.c $(DEFV5) $(INCGRO) -I$(SVM) -I$(GKUT)/include

gkut:
	make CC=$(CC) CFLAGS=$(MCFLAGS) GROMACS=$(GROMACS) VGRO=$(VGRO) -C $(GKUT)

clean:
	make clean -C $(SVM) && make clean -C $(GKUT) \
	&& rm -f $(BUILD)/*.o $(BUILD)/g_ensemble_res_comp

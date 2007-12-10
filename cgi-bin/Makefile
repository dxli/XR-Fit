PHONY: all

GXX=g++
LDFLAGS=-lpthread -lmysqlpp -lgrace_np -lgsl -lgslcblas -lgrace_np -lgomp -Wl,-O1 -Wl,-as-needed
INCLUDE_DIRS=-I/usr/include/mysql++ -I/usr/include/mysql
CXXFLAGS=-O3 -march=nocona -m64 -fexpensive-optimizations -fopenmp
SRCS=GARealGenome.cpp multilayer.cpp logUtility.cpp Task.cpp nph_mysql.cpp elements.cpp file_copy.cpp
gixosSRCS=GARealGenome.cpp gixos_multilayer.cpp logUtility.cpp Task.cpp gixos_nph_mysql.cpp elements.cpp file_copy.cpp gixos.cpp
SRCS2=Task.cpp nph_mysql.cpp elements.cpp logUtility.cpp

all: nph dispatcher plotmc1 plot_gixos gixos_nph

nph: ${SRCS} Makefile nph.cpp defines.h
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SRCS} nph.cpp

dispatcher: ${SRCS2} Makefile dispatcher.cpp
	${GXX} ${LDFLAGS} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ ${SRCS2} dispatcher.cpp

plotmc1: plotmc1.cpp Makefile
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ $< ${LDFLAGS}

gixos_nph: ${gixosSRCS} Makefile gixos_nph.cpp defines.h
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} ${LDFLAGS} -o $@ ${gixosSRCS} gixos_nph.cpp

plot_gixos: plot_gixos.cpp Makefile
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ $< ${LDFLAGS}

clean:
	rm *.o nph gixos_nph plotmc1 dispatcher nph.log -f


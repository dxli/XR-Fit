PHONY: all

<<<<<<< HEAD:cgi-bin/Makefile
GXX=g++
LDFLAGS=-lpthread -lmysqlpp -lgrace_np -lgsl -lgslcblas -lgrace_np -lgomp -Wl,-O1 -Wl,-as-needed
=======
CXX =		g++
CXXFLAGS =	-g -O2  -DHAVE_CONFIG_H
LDFLAGS =	 -lpthread -lmysqlpp -lgslcblas -lgsl -lgrace_np -lgomp 
#LDFLAGS=-lpthread -lmysqlpp -lgrace_np -lgsl -lgslcblas -lgrace_np -lgomp -Wl,-O1 -Wl,-as-needed
>>>>>>> origin/origin/master:cgi-bin/Makefile
INCLUDE_DIRS=-I/usr/include/mysql++ -I/usr/include/mysql
<<<<<<< HEAD:cgi-bin/Makefile
CXXFLAGS=-O3 -march=nocona -m64 -fexpensive-optimizations -fopenmp
=======
#CXXFLAGS=-O3 -march=nocona -m64 -fexpensive-optimizations -fopenmp
>>>>>>> origin/origin/master:cgi-bin/Makefile
SRCS=GARealGenome.cpp multilayer.cpp logUtility.cpp Task.cpp nph_mysql.cpp elements.cpp file_copy.cpp
gixosSRCS=GARealGenome.cpp gixos_multilayer.cpp logUtility.cpp Task.cpp gixos_nph_mysql.cpp elements.cpp file_copy.cpp gixos.cpp
SRCS2=Task.cpp nph_mysql.cpp elements.cpp logUtility.cpp

all: nph dispatcher plotmc1 plot_gixos gixos_nph

nph: ${SRCS} Makefile nph.cpp defines.h
<<<<<<< HEAD:cgi-bin/Makefile
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SRCS} nph.cpp
=======
	${CXX} ${INCLUDE_DIRS} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SRCS} nph.cpp
>>>>>>> origin/origin/master:cgi-bin/Makefile

dispatcher: ${SRCS2} Makefile dispatcher.cpp
<<<<<<< HEAD:cgi-bin/Makefile
	${GXX} ${LDFLAGS} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ ${SRCS2} dispatcher.cpp
=======
	${CXX} ${LDFLAGS} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ ${SRCS2} dispatcher.cpp
>>>>>>> origin/origin/master:cgi-bin/Makefile

plotmc1: plotmc1.cpp Makefile
<<<<<<< HEAD:cgi-bin/Makefile
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ $< ${LDFLAGS}
=======
	${CXX} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ $< ${LDFLAGS}
>>>>>>> origin/origin/master:cgi-bin/Makefile

gixos_nph: ${gixosSRCS} Makefile gixos_nph.cpp defines.h
<<<<<<< HEAD:cgi-bin/Makefile
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} ${LDFLAGS} -o $@ ${gixosSRCS} gixos_nph.cpp
=======
	${CXX} ${INCLUDE_DIRS} ${CXXFLAGS} ${LDFLAGS} -o $@ ${gixosSRCS} gixos_nph.cpp
>>>>>>> origin/origin/master:cgi-bin/Makefile

plot_gixos: plot_gixos.cpp Makefile
<<<<<<< HEAD:cgi-bin/Makefile
	${GXX} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ $< ${LDFLAGS}
=======
	${CXX} ${INCLUDE_DIRS} ${CXXFLAGS} -o $@ $< ${LDFLAGS}
>>>>>>> origin/origin/master:cgi-bin/Makefile

clean:
	rm *.o nph gixos_nph plotmc1 dispatcher nph.log -f


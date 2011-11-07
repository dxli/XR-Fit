#general c++ compiling rule
bindir = /var/www/localhost/cgi-bin
AM_CPPFLAGS = -fomit-frame-pointer -fexpensive-optimizations -DMYSQLPP_MYSQL_HEADERS_BURIED -lomp
bin_PROGRAMS =  nph dispatcher plotmc1 plotgixos gixosnph
nph_SOURCES = nph.cpp
nph_LDADD = GARealGenome.o multilayer.o logUtility.o Task.o nph_mysql.o elements.o file_copy.o
dispatcher_SOURCES = dispatcher.cpp
dispatcher_LDADD = nph_mysql.o logUtility.o elements.o Task.o
gixosnph_SOURCES = gixosnph.cpp
gixosnph_LDADD = GARealGenome.o gixosmultilayer.o logUtility.o Task.o gixosnph_mysql.o elements.o file_copy.o gixos.o 
plotmc1_SOURCES = plotmc1.cpp
plotgixos_SOURCES = plotgixos.cpp
INCLUDES = -I/usr/include/mysql++

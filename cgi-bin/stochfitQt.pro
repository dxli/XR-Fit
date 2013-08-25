#-------------------------------------------------
#
# Project created by QtCreator 2013-08-25T02:16:16
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = stochfitQt
TEMPLATE = app

LIBS += -lgsl -lgslcblas -lm

SOURCES += main.cpp\
        mainwindow.cpp \
    multilayer.cpp \
    GARealGenome.cpp

HEADERS  += mainwindow.h \
    multilayer.h \
    GARealGenome.h

FORMS    += mainwindow.ui

#pragma once



#include <string>
#include <cstring>
#include <mysql++.h>
#include <iostream>
#include <cstdlib>
#include "gixosnph.h"
#include "Task.h"

using namespace std;
using namespace mysqlpp;

#define DATABASE "nph_data"
#define HOST     "localhost"
#define USER     "nph"
#define PASSWORD "nph"


void getTask(Task * t, int id);
int getNextTaskId();
int getRunning(int);
int getStarted(int);
int updateProgress(int id, int progress,double);
int updateStarted(int id, int started);
int updateDone(int id, bool done);
int updateRunning(int id, int running);
int updateError(int id, int error);


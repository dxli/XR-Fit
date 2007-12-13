#pragma once

/* ----------------------------------------------------------------------------
// Modified to do reflectivity fitting
// au.gafitm.cpp should be optimized to do monolayer fitting
// restrict the beta profile by relating it to the density profile

rho= k ( rho_a (1-x) + rho_b x)
beta= k ( beta_a (1-x) + beta_b x)

modified to do nph-cgi

---------------------------------------------------------------------------- */
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <complex>
#include <cctype>

#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
 
#include "GARealGenome.h"
#include "gixosmultilayer.h"
#include "defines.h"
#include "logUtility.h"

#include "Task.h"
#define N_Av 6.02214179e23

using namespace std;

ostream& operator << (ostream& os, GARealGenome genome);
void read_in_data(string * data0, char ** argv, int argc);
void getData(int id, string *data0, Task *t);
void run(int taskId);
int file_copy(string,string);


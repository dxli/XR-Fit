#ifndef GARealGenome_H
#define GARealGenome_H

#include <vector>
#include <cmath>
#include <pthread.h>


#include "defines.h"

using namespace std;

class GARealGenome
{
private:
    vector <double> genome;
    double dgene2;

public:

    double gene_min,gene_max;
    double dgene;
    double y_factor;

    GARealGenome(int l,double min,double max);
    void copy(GARealGenome);
    int size();
    double gene(int i);
    void gene(int i,double x);
    void adjust(int,double);

    //single random jumping of the genome
    void randomize0(double x);

    //random jumping of the genome
    void randomize(double x);
};
#endif


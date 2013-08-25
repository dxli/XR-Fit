#include "GARealGenome.h"
#include "logUtility.h"
#include <cstdlib>


GARealGenome::GARealGenome(int l,double min,double max)
{
	genome.resize(l);
	gene_min=min;
	gene_max=max;
	dgene=max-min;
	dgene2=dgene+dgene;
	//cout<<"gene_min= "<<gene_min<<", gene_max="<<gene_max<<endl;
}

void GARealGenome::copy(GARealGenome b){
this->y_factor=b.y_factor;
this->genome.resize(b.size());
vector<double>::iterator pgene=this->genome.begin();
vector<double>::iterator pgeneb=b.genome.begin();
while(pgene != this->genome.end()) {
	*pgene++ =fabs(fmod(*pgeneb++ +dgene,dgene2)-dgene);
	}
	/*
this->gene_min=b.gene_min;
this->gene_max=b.gene_max;
this->dgene=b.dgene;
*/
}

int GARealGenome::size()
{
	return(genome.size());
}


double GARealGenome::gene(int i)
{
	return genome.at(i);
}

void GARealGenome::gene(int i,double x)
{
	genome.at(i)=x;
}


void GARealGenome::adjust(int jj, double x)
{ //single random jumping of the genome
	genome.at(jj) =fabs(fmod(genome.at(jj)+x +dgene,dgene2)-dgene);
	}

void GARealGenome::randomize(double /*x*/)
{ //random jumping of the genome
	//genome.at(0) *= exp((random()/(RAND_MAX+1.0)-0.5)*x);
	double a=0.15*dgene;
    for(unsigned jj=0;jj<genome.size();jj++) {
        genome.at(jj) =fabs(fmod( genome.at(jj) +a*(random()/(RAND_MAX+1.0)-0.5) +dgene,dgene2)-dgene);
	}
};


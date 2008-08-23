#pragma once
#include <vector>
#include <complex>
#include <fstream>
#include <istream>
#include <iterator>
#include <pthread.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <omp.h>


#include "GARealGenome.h"
#include "defines.h"
#include "logUtility.h"
#include "xrdata.h"


using namespace std;

class multilayer
{
private:
    vector <complex<double> >  nk,dzk;
    vector <double> rhok,dz;
    vector <double> genome0,genome1;
    //mkdensity parameter
    //double gene_ra,gene_rb,gene_rc,gene_ba,gene_bb,gene_bc;
public:
    double dz0,dz1,slab,sdyi;
    double alph0,dth,sigma0,gamma;
    //vector <double> xi,yi,dyi,thetai,x2i;
    vector <xrdata> ref0,ref1;
    string fnpop,fnref,fnrf,fnrho;
    unsigned int nl,glength,dside;
    //double rho_a,rho_b,beta_a,beta_b;
    complex<double> n_bulk0,n_bulk1,n_bulk12;
    //double rho0,beta0,rho12,beta12;
    double lambda,k0;
    complex<double> kdz0; //k0*dz0*2i
    double qmin,qmax;//qz range for fitting
    multilayer(unsigned int ,double ,double,double);
    multilayer(){}
    // ~multilayer();
    void readref(string);
    double rf(vector<xrdata>::iterator);
    double fresnel(double t);
    complex<double> rho_gene(double t);
    double objective(GARealGenome  *);
    void genomerf(GARealGenome *);
    void genomerf(GARealGenome *,double,double,int);
    void mkdensity(GARealGenome *);
    double cT(double),yet_qz(double),sip_qz(double);
    void setbulk(double ,double ,double ,double ,double ,double );
};


/* ----------------------------------------------------------------------------
  mbwall 5jan96
  Copyright (c) 1995-1996  Massachusetts Institute of Technology
// Modified to do reflectivity fitting
// au.gafitm.cpp should be optimized to do monolayer fitting
// restrict the beta profile by relating it to the density profile

rho= k ( rho_a (1-x) + rho_b x)
beta= k ( beta_a (1-x) + beta_b x)
modified to do si-pb interface
---------------------------------------------------------------------------- */
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <istream>
#include<vector>
#include<sstream>
#include<string>
#include<algorithm>
#include<iterator>

#define rhoE_maximum 5.
#define slab_length 55.
#define MC_STEPSIZE 0.25
#define sg_factor       100.
#define intp_factor     4
#define density_factor  10.
#define intg_factor     256.

// Multilayer reflection and transmission
// Parrett scheme, L. G. Parrett, Phys. Rev. 95(2), 359(1954)
// formulaa to use
//      n \cos(\theta) = \cos(\theta_0)
//
//      (E_{k,+} + E_{k,-} ) \cos(\theta_k) = ( E_{k+1,+}/a + E_{k+1,-} a) \cos( \theta_{k+1})
//      (E_{k,+} - E_{k,-} ) \sin(\theta_k) = ( E_{k+1,+}/a - E_{k+1,-} a) \sin( \theta_{k+1})
//      a = exp( i k_{k+1} d_{k+1})
//      k_{k+1} = sqrt( n^2 - \cos^2(\theta_0)) k_0

using namespace std;
#include<complex>
#include<vector>

class GARealGenome {
private:
    vector <double> genome;
    double dgene;
public:
    double gene_min,gene_max;
    int length;
    GARealGenome(int l,double min,double max){
        genome.resize(l);
        length=l;
        gene_min=min;
        gene_max=max;
        dgene=max-min;
        //cout<<"gene_min= "<<gene_min<<", gene_max="<<gene_max<<endl;
    };
    int size(){
        return(length);
    };
    double gene(int i){
        return genome.at(i);
    };
    void gene(int i,double x){
        genome.at(i)=x;
    };
    void randomize0(double x){ //single random jumping of the genome
        int jj=(int) ( (double) random()/(1.0+RAND_MAX)*genome.size());
        if(jj) {
                genome.at(jj) =gene_min+(genome.at(jj)-gene_min)*exp((random()/(RAND_MAX+1.0)-0.5)*x);
            if(genome.at(jj)<gene_min) {
                genome.at(jj)=gene_min+gene_min-genome.at(jj);
                return;
            }
            genome.at(jj) *=exp((random()/(RAND_MAX+1.0)-0.5)*x);
            if(genome.at(jj)>gene_max) {
                genome.at(jj)=gene_max+gene_max-genome.at(jj);
            }
        }
    };
    void randomize(double x){ //random jumping of the genome
        for(int jj=1;jj<genome.size();jj++) {
            double a= gene_min+(genome.at(jj)-gene_min)*exp((random()/(RAND_MAX+1.0)-0.5)*x);
            if(a<gene_min) genome.at(jj)= gene_min+gene_min -a;
            else if(a>gene_max) genome.at(jj)= gene_max+gene_max -a;
                else genome.at(jj)=a;
        }
        genome.at(0) *= exp((random()/(RAND_MAX+1.0)-0.5)*x);
    };
};

ostream& operator << (ostream& os, GARealGenome genome) {
    // Output the atom position to the stream in the format n/d
    for(int jj=0;jj<genome.size();jj++) os<<' '<<genome.gene(jj);
    return os;
}


class multilayer {
private:
    vector <complex<double> >  thetak,kk,ak,nk,ep,em;
    vector <double> rhok,dz;
    vector <double> genome0,genome1;
public:
    double dz0,dz1,sdyi;
    double smooth_factor;
    vector <double> xi,yi,dyi,thetai,x2i;
    string fnpop,fnrf,fnrho;
    unsigned int nl,glength,side;
    double rho_a,rho_b,beta_a,beta_b;
    complex<double> n_bulk0;
    double lambda,k0;
    multilayer(unsigned int i,double l,vector<double> d0,vector <double> rho0,vector <double> beta);
    multilayer(){
    }
    // ~multilayer();
    void init(unsigned int i,double l,vector<double> d0,vector <double> rho0,vector <double> beta);
    void readref(string fn);
    void theta(double t);
    double rf(double t);
    double rho_conv(double x);
    double objective(GARealGenome  *);
    void genomerf(GARealGenome *);
    void mkdensity(GARealGenome *);
    void setbulk0(double a,double b){
        n_bulk0= complex<double>(2.8179e-5*a*lambda*lambda/(2.*M_PI),-b);
    };
    void setbulk1(double a, double b){
        nk[nl-1]=
            complex<double>(2.8179e-5*a*lambda*lambda/(2.*M_PI),-b);
    };
};

multilayer::multilayer( unsigned int i, double l,vector<double> d0,vector<double> rho0,vector<double> beta)
{
    nl= i;
    lambda=l;
    k0=2*M_PI/lambda;

    if(rho0.size() < i) {
        //cout<<"Not enough rho0 data"<<endl;
        exit(0);
    }
    //rhok.resize(nl);
    //dz.resize(nl);
    rhok = rho0;
    dz = d0;
    dz0=d0[1];
    thetak.resize(nl);
    kk.resize(nl);
    nk.resize(nl);
    ak.resize(nl);
    ep.resize(nl);
    em.resize(nl);
    int j;
    for(j=0;j<nl;j++){
        nk[j]= complex<double> (2.8179e-5*rho0[j]*lambda*lambda/(2.*M_PI),0.-beta[j]);
    }
    nk[0]=0.;
}

void multilayer::init( unsigned int i, double l,vector<double> d0,vector<double> rho0,vector<double> beta)
{
    nl= i;
    lambda=l;
    k0=2*M_PI/lambda;

    if(rho0.size() < i) {
        //cout<<"Not enough rho0 data"<<endl;
        exit(0);
    }
    //rhok.resize(nl);
    //dz.resize(nl);
    rhok = rho0;
    dz = d0;
    dz0=d0[1];
    thetak.resize(nl);
    kk.resize(nl);
    nk.resize(nl);
    ak.resize(nl);
    ep.resize(nl);
    em.resize(nl);
    for(unsigned j=0;j<nl;j++){
        nk[j]= complex<double> (2.8179e-5*rho0[j]*lambda*lambda/(2.*M_PI),0.-beta[j]);
    }
    nk[0]=0.;
}

void multilayer::readref( string fn)
{
    double x,y,dy;
    sdyi=0.;

        istringstream iss (fn);
        vector < double >va;
        std::copy (istream_iterator < double >(iss),
                   istream_iterator < double >(), back_inserter (va));
        int jj= (va.size()/3)*3;
    for(int ii=0;ii<jj;ii+=3)
    {
        x=va.at(ii);
        y=va.at(ii+1);
        dy=va.at(ii+2);

        if(fabs(y) > 1.e-18){ // we drop point zero, since we use relative errorbars
            //                        dy=sqrt(dy);
            //if(x<0.15) dy *= 0.1;
            //disable weights
            //dy=y;
            xi.push_back(x);
            yi.push_back(y);
            //cout<<x<<' '<<dy/y<<endl;
            double wt=2./(1.+abs(x/0.2))/(0.05+dy/y);
            /*
            wt=y/dy;
            if(xi.size()>1) {
                wt= 1./(1./wt + 10./sqrt(dyi[0]));
            }
            */
            wt *= wt;
            dyi.push_back(wt);
            sdyi += wt;
            x2i.push_back(-x*x*sg_factor);
            thetai.push_back(asin(x*lambda/(4*M_PI)));
        }
    }
    sdyi = 1./sdyi;
    //cout<<xi.size()<<"\t"<<sdyi<<endl;
}

void multilayer::theta( double t)
{
    thetak[0]=t;
    unsigned int j;
    complex<double> a0,b0; // (1 +nk[0]) *(1- t^/2+t^4/24) = 1 -(t*t/2 +nk[0]*t*t/2-nk[0])
    //a0 =  t*t/2. - t*t*t*t/24. + nk[0]*t*t/2. - nk[0];
    double c=t*t;
    a0 =  c*(1./2.-c*(1./24.-1./720.*c)); // 1 - cos(t)
    for(j=1;j<nl;j++){
        // (1 - a0) / (1 - nk[j]) = 1 - ( a0 -nk[j] +a0*nk[j]-nk[j]*nk[j])
        b0=a0 - nk[j] +a0 * nk[j]+ nk[j]*nk[j];
        complex<double> c0;
        c0=sqrt(b0+b0); // to make b0^2/2
        thetak[j]= c0*(1. +b0*( 1./12.+b0*(3./160.+5./896.*b0))); // acos(1 - c0^2/2)
    }
    complex<double> i2pi(0.,1.);
    for(j=1;j<nl-1;j++){
        //[ (1-nk[j]) ]^2 - cos(t)^2 = t*t/2 - 2*nk[j] + nk[j]^2
        //kk[j]=k0*sqrt( t*t/2. -t*t*t*t/24. -2.*nk[j] +nk[j]*nk[j]);
        kk[j]=k0*(1. - nk[j])*sin(thetak[j]);
        ak[j]=exp( i2pi* kk[j]*dz0);
    }
    kk[0]=k0;
    ak[0]=ak[j]=1.;
}

double multilayer::rf(double t)
{
    int i;
    theta(t);
    ep[nl-1]=1.;
    em[nl-1]=0.;
    for(i=nl-2;i>=0;i--){
        complex<double> a,ep0,em0,c,d;
        a=ak[i+1]; // ak[nl-1]=1.
        ep0=ep[i+1]/a;em0=em[i+1]*a;
        c=(ep0+em0)*cos(thetak[i+1])/cos(thetak[i]);
        d=(ep0-em0)*sin(thetak[i+1])/sin(thetak[i]);
        ep[i]=0.5*(c+d);
        em[i]=0.5*(c-d);
    }
    double a=abs(em[0]/ep[0]);

    return( a*a);
}

void multilayer::genomerf(GARealGenome * g)
// output for plotting
{
    //if(g0.length() < nl+nl -3) {
    //cout<<"Incorrect genome size "<<g0.size()<<endl;
    //}
    unsigned int j,jj=nl-1;
    //cout<<dz0<<endl;
    ofstream out(fnrho.c_str());
    mkdensity(g);
    double x=dz0;
    for(j=0;j<=jj;j++) {
        double r0;
        int j0= (int) ( x/dz1+0.5-side)+1;
        if(j0<1) r0=0.;
        else if (j0<glength) r0=g->gene(j0);
        else r0=1.0;
        r0=rho_a - r0*rho_b;
        out<<x<<' '<<nk[j].real()/nk[jj].real()<<' '<<r0<<' '<<nk[j].imag()/nk[jj].imag()<<endl;
        x+=dz0;
    }
    out.close();
    out.open(fnrf.c_str());

    double x0=xi[0];
    double x1=xi[xi.size()-1];
    double dx=(x1-x0)/300.;
    x1=1.1*x1-0.1*x0;
    for(x=x0;x<=x1;x += dx){
        // cout<<b*180./M_PI<<" "<<a<<endl;
        // weighted using dyi[j], E(y) =  yi / dyi,
        out<<x<<" "<<rf(asin(x*lambda/(4*M_PI)))<<endl;
    }
    out.close();
}

void multilayer::mkdensity(GARealGenome * g)
{
    unsigned int j,jj=nl-1;
    double b=0.,sy3=0.;
    double x=dz0,s0=sg_factor*g->gene(0),r0=0.,i0=0.;
    double s1=-0.5/s0,ds=sqrt(s0),dx=ds/intg_factor;
    double s2=sqrt(1./(2*M_PI))/intg_factor,x0, x1= 4.*ds;
    unsigned int jl,jr;
    x=(side*dz1 -x1)/dz0;
    if(x>0) jl= (unsigned int) ( x+0.5); else jl=1;
    jr= (unsigned int) ( (side+glength+1)*dz1/dz0 +0.5);
    if (jr<jj) jr=jj;
    //   cout<<"dz0="<<dz0<<" dz1="<<dz1<<" side="<<side<<" glength="<<glength<<endl;

    //  cout<<"j1="<<jr<<endl;

    for(j=0;j<jl;j++){
        nk[j]= n_bulk0;
    }
    x= j*dz0;
    for(;j<jr;j++){
        //convolution
        r0=i0=0.;
        x0= -x1;
        while(x0<x1){
            double a = exp(s1*x0*x0);
            int j0= (int) ( (x+x0)/dz1+1.5-side);
            if (j0 >=1) {
                if(j0>glength) {
                    r0 += a;
                    i0 += a;
                }else{
                    r0 += a*(rho_a - rho_b*g->gene(j0));
                    i0 += a*(beta_a - beta_b*g->gene(j0));
                }
            } else {
                r0 += a*rho_a;
                i0 += a*beta_a;
            }
            x0+=dx;
        }
        r0 *=s2;
        i0 *=s2;
        nk[j]= complex<double> (r0*nk[jj].real(),i0*nk[jj].imag());
        x+=dz0;
    }
    for(;j<jj;j++){
        nk[j]= nk[jj];
    }


}

double multilayer::objective(GARealGenome * g)
// objective function
{
    //if(g0.length() < nl+nl -3) {
    //cout<<"Incorrect genome size "<<g0.size()<<endl;
    //}
    //        for(j=0;j<nl;j++) cout<<j<<" "<<nk[j]<<endl;
    mkdensity(g);
    double sy=0.,sy2=0.,a;
    //double d= g0.gene(0);
    //for(unsigned int j=0;j<xi.size() && xi[j] <=smooth_factor;j++){
    for(unsigned int j=0;j<xi.size();j++){
        //a = exp(d*x2i[j])* ( 1./6.*(rf(b - 0.005*M_PI/180.)+rf(b + 0.005*M_PI/180.))+4./6.*rf(b));
        double a = rf(thetai[j]);
        //cout<<thetai[j]<<" , "<<a<<endl;
        // weighted using dyi[j], E(y) =  yi / dyi,
        //               disabling weight
        a= log(a/yi[j]);
        //sy += a*dyi[j];
        sy2 += a*a*dyi[j];
    }

    //sy *= sdyi;
    sy2 *= sdyi;
    //sy2 -= sy*sy;
    //a= sy*sy*(sy2 - sy*sy);
    //cout<<"ob: "<<a<<endl;
    return((double) sy2);
}

multilayer ml0;
int
main(int argc, char** argv)
{
    int i=2048; // number of layers in Parrett reflectivity
    // e=1.602176e-19
    // h=6.626069e-34
    // c=2.997925e8
    double l0,e_x=7.1e4 /*in eV now*/;
    {
        double  e=1.602176e-19, h=6.626069e-34,  c=2.997925e8;
        l0= h*c/(e_x*e)*1.e10;
    }
    //cout<<"x-ray LAMBDA= "<<l0<<" angstrom\n";

    double Nav=6.022e23;
    vector<double> rho0,beta,dz;
    double rhoBulk1=10.66,rhoBulk0=2.33; // Gold density 6.022 10^23 19.3 1000. / 0.1974 79  10^-30;
    double a_Bulk1=207.2,a_Bulk0=28.09;
    double z_Bulk1=82,z_Bulk0=14;
    rhoBulk1 *= Nav*1e6/a_Bulk1*z_Bulk1 *1.e-30;
    rhoBulk0 *= Nav*1e6/a_Bulk0*z_Bulk0 *1.e-30;
    double betaBulk1=5.13e-9*10.66/11.35,betaBulk0=8.033e-11,beta01=2*betaBulk1-betaBulk0;
    rho0.resize(i);
    beta.resize(i);
    dz.resize(i);
    unsigned int j;
    double dz0=2.5*slab_length/(i-2);
rhoBulk0=0.;
betaBulk0=0.;
    ml0.rho_a=rhoBulk0/rhoBulk1;
    ml0.rho_b=ml0.rho_a -1.;
    ml0.beta_a=betaBulk0/betaBulk1;
    ml0.beta_b=ml0.beta_a -1.;
    for(j=1;j<i-1;j++) {
        dz[j]= dz0;  // i slabs, i-1 interfaces, i - 2 divided
        rho0[j]=rhoBulk0;
        beta[j]=betaBulk0;
    }
    rho0[0]=0.;
    beta[0]=0.;
    rho0[i-1]=rhoBulk0;
    beta[i-1]=betaBulk0;
    dz[0]=dz[i-1]=0.;
    //cout<<"i="<<i<<" l="<<l0<<endl;
    ml0.init(i,l0,dz,rho0,beta);
    ml0.setbulk0(rhoBulk0,betaBulk0);
    ml0.setbulk1(rhoBulk1,betaBulk1);
    ml0.smooth_factor=0.05;
    vector<double> g0;
    j=(unsigned int ) (3.*4./( pow((double) 2,intp_factor)*ml0.dz0)+0.5);
    ml0.side=j;
    unsigned int glength= (i>>intp_factor) - 2 - j-j ;
    ml0.glength=glength;
    g0.resize(glength +1);
    for(j=1;j<=glength;j++){
        g0[j]=1.;
    }
    g0[0]=2.;
    g0[0] = g0[0]* g0[0]/(sg_factor);
        cout<<"Content-Type:text/html;charset=iso-8859-1\r\n\r\n";
        cout<<"<html><head><title>Results</title></head><body>\n";
string data0;
        {
                char *lenstr;
                lenstr = getenv("CONTENT_LENGTH");
                if(lenstr==NULL) exit(0);
int data_len;
        if( sscanf(lenstr,"%d",&data_len)!=1) exit(0);
        if( data_len<=0 || data_len>=10000) exit(0);
int ii=0;
while(ii++<data_len) {
        char c0=getchar();
        if(c0 == '=') break;
}
while(ii++<data_len) {
        char c0=getchar();
        if(c0=='%') {
                ii+=2;
                getchar();getchar();
                c0='\n';
                data0 += "<br>";
        }
        if(c0=='+') c0=' ';
        data0.push_back(c0);
}
        }
        cout<<data0<<"<br></body></html>";
    ml0.readref(data0);
    string fnref("tmpout");
    ml0.fnpop=fnref+string("pop.dat");
    ml0.fnrf=fnref+string("rf.dat");
    ml0.fnrho=fnref+string("rho.dat");
    ml0.dz1 =(( (unsigned int) 1)<<intp_factor)*ml0.dz0;
    cout<<"dz0="<<ml0.dz0<<" dz1="<<ml0.dz1<<endl;

    //cout << "Fitting using galib, http://lancet.mit.edu/ga/\n\n";
    cout.flush();

    // See if we've been given a seed to use (for testing purposes).  When you
    // specify a random seed, the evolution will be exactly the same each time
    // you use that seed number.
    char random_buf[257];
    {
        string sr("/dev/urandom");
        struct timeval time0;
        gettimeofday(&time0,NULL);
        ifstream in1(sr.c_str());
        if(! in1.is_open()) {
            cout<<"Can not open "<<sr<<endl;
            srandom(time0.tv_usec);
        }else {
            in1.read(random_buf,256*sizeof(char));
            in1.close();
            initstate(time0.tv_usec,random_buf,256);
            setstate(random_buf);
        }
    }
    //cout<<"random: "<<random()<<endl;


    ofstream outfile;

    //  GARealAlleleSet alleles(0., MAX_VALUE);
    double agene_min=ml0.rho_a/ml0.rho_b;
    double agene_max=( ml0.rho_a -rhoE_maximum)/ml0.rho_b;
    cout<<"gene: [ "<<agene_min<<","<<agene_max<<"]\n";
    GARealGenome genome(1+glength , agene_min,agene_max);
    genome.randomize(MC_STEPSIZE);
    /*
        for(int jj=0;jj<genome.size();jj++) genome.gene(jj,(double) (g0[jj]* exp( -0.01*random()/(RAND_MAX+1.0))));
        */
    // dump the initial population to file

    //cout << "printing initial population to file..." << endl;

    GARealGenome genome1=genome;

    ifstream infile(ml0.fnpop.c_str());

    int kk=1;
    if(infile.is_open()){
        double xi,yi,iyi;
        infile>>xi>>yi>>iyi;
        int jj=0;
        genome1.gene(jj++,iyi/sg_factor);

        while(jj< genome1.size()){
            if(infile>>xi>>yi>>iyi){
                genome1.gene(jj,yi);
                jj++;
            } else {
                kk=0;
                break;
            }
        }
        infile.close();
    } else kk=0;
    if(kk) {
        cout<<"Using the saved genome in "+ml0.fnpop<<endl;
    }else{
        cout<<"Using a random genome"<<endl;
        for(int jj=1;jj<genome1.size();jj++)
            genome1.gene(jj,agene_min+ (double) random()/(1.0+RAND_MAX)*(agene_max-agene_min));
        genome1.gene(0,random()/(1.0+RAND_MAX)*0.15+0.05);
    }


    genome=genome1;
    {
        double agene1=genome.gene(0);
        if(agene1>0.2) agene1=0.2;
        agene1=fabs(agene1);
        genome.gene(0,agene1);
    }

    cout<<genome<<endl;
    //ml0.genomerf(genome);
    //return 0 ;

    int ii=0;
    double score0,score1,score2,score20;

    GARealGenome genome2=genome;
    score20=score2=score0=ml0.objective(&genome);
    double mc_step=MC_STEPSIZE;
    unsigned acc0=0;
    unsigned int jmax=genome.size()*4;
    for(unsigned isteps=1;;isteps++){
        for(unsigned int jj=0; jj<jmax;jj++){
            genome1=genome;
            genome1.randomize0(mc_step);
            score1=ml0.objective(&genome1);
            if( score1< score0) {
                acc0++;
                genome=genome1;
                score0=score1;
            }
        }
        double acpt=(double) acc0/jmax;
        acc0=0;
        if(score0<score20){
            genome2=genome;
            score20=score0;
        }
        cout<<isteps<<" acpt="<<acpt<<' '<<score0<<'('<<score2<<","<<score20<<") "<<mc_step<<endl;
        if(acpt>0.3) mc_step *=1.2;
        if( acpt< 0.2 || score0 == score2) {
            //if( ml0.smooth_factor <0.55) ml0.smooth_factor +=0.02;
            if(mc_step<5e-3 && fabs(score0 - score2) <=1e-5) {
                mc_step=MC_STEPSIZE;
               // cout<<genome<<endl;
                genome.randomize(mc_step);
                score0=ml0.objective(&genome);
                //cout<<"after randomize():\n"<<genome<<endl;
            } else mc_step*=0.8;
        } else score2=score0;
        outfile.open(ml0.fnpop.c_str());
        int jj=0;
        double x=0.,dx=(double) ml0.dz1;
        outfile<<ii<<' '<<score20<<' '<<sg_factor*genome2.gene(jj++)<<endl;
        do{
            outfile<<x<<' '<<genome2.gene(jj)<<" 0.\n";
            x += dx;
            jj++;
        }while(jj<genome2.size());

        outfile.close();
        ml0.genomerf(&genome2);
    }
    return 0;
}




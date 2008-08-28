#include "multilayer.h"
#include <iostream>
#include "nph.h"
using namespace std;

multilayer::multilayer( unsigned int i, double l,double dz00, double dz10)
{
    nl= i;
    lambda=l;
    k0=2*M_PI/lambda;
    nk.resize(nl);
    dz0=dz00;
    slab=(nl-2)*dz00;
    dz1=dz10;
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
    ofstream out1(fnref.c_str());
    int jj= (va.size()/3)*3;
    for (int ii=0;ii<jj;ii+=3)
    {
        x=fabs(va.at(ii));
        y=va.at(ii+1);
        dy=va.at(ii+2);
        /*
        wt=y/dy;
        if(xi.size()>1) {
            wt= 1./(1./wt + 10./sqrt(dyi[0]));
        }
        */
        double sink0=x*lambda/(4*M_PI);
        double fresnel0=1./fresnel(sink0);

        out1<<x<<' '<<y*fresnel0<<' '<<dy*fresnel0<<endl;

        if (y > 1.e-18 && x>=qmin && x<=qmax) { // we drop point zero, since we use relative errorbars
            double wt=1./(0.075+dy/y);
            wt *= wt;
            sdyi += wt;

            //                        dy=sqrt(dy);
            //if(x<0.15) dy *= 0.1;
            //disable weights
            //dy=y;
            //cout<<x<<' '<<dy/y<<endl;
            //double wt=2./(1.+abs(x/0.05))/(0.3+dy/y);
            //cout<<x<<" "<<fresnel0<<endl;
            //a0 =  t*t/2. - t*t*t*t/24. + nk[0]*t*t/2. - nk[0];
            //double c=t*t;
            //double cosk0=  c*(1./2.-c*(1./24.-1./720.*c)); // 1 - cos(t)
            complex<double> a0=n_bulk0 + sink0*sink0;
            //a0 += n_bulk0*(1.0 - a0);
            ref0.push_back(xrdata(x,y*fresnel0,dy*fresnel0,y,wt,sink0,fresnel0,sink0,a0));
        }
    }
    out1.close();
    sdyi = 1./sdyi;
    //cout<<xi.size()<<"\t"<<sdyi<<endl;
    //write out reflectivity data
    // cout<<"Writing: "<<fnref<<endl;
    string::size_type locr=fnref.rfind(".",fnref.size());
    string::size_type locl=fnref.rfind("-",fnref.size());
    if (locr != string::npos) {
        string fnref0=string("tmpout-")+fnref.substr(locl+1, locr-locl-1)+string("-ref.dat");
        cout<<"cp "<<fnref<<" "<<fnref0<<endl;
        if (file_copy(fnref, fnref0)) {
            cout<<"Error: cp "<<fnref<<" "<<fnref0<<endl;
        }
    }
    //ak.at(0)=1.;
    kdz0=complex<double>(0,2)*k0*dz0;
    //xi for plotting rf()
    //double x0=fabs(ref0.at(0).xi)*0.2;
    double x0=0.01;
    double x1=ref0.at(ref0.size()-1).xi;
    double dx=(x1-x0)/300.;
    x1=1.1*x1-0.1*x0;
    for (x=x0;x<=x1;x += dx) {
        // cout<<b*180./M_PI<<" "<<a<<endl;
        // weighted using dyi[j], E(y) =  yi / dyi,
        double t=x*lambda/(4*M_PI);
        complex<double> a0=n_bulk0+t*t;
        double fresnel0=1./fresnel(t);
        ref1.push_back(xrdata(x,1.,1.,1.,1.,t,fresnel0,t,a0));
    }

}



double multilayer::fresnel(double t)
//fresnel reflectivity
{
//    double cosk0=t*t;
//   cosk0 =  cosk0*(1./2.-cosk0*(1./24.-1./720.*cosk0)); // 1 - cos(t)
    //double cosk0=1-a;
    //c=sqrt(a+a);
    //double sink0=sqrt(cosk0*(2.0-cosk0));
    complex<double> sink0 =n_bulk0+t*t;
    complex<double> sink1 =sink0 - n_bulk1;
    sink0=sqrt(sink0);
    sink1=sqrt(sink1);
    //cosk0=1-cosk0;
    //sink1=(sink1 - n_bulk1)*(1.+n_bulk1); // (1 -a0) /(1-nk)=(1-a0)(1+nk+nk^2)=1 -(a0 +a0nk-nk-nk^2)
    //complex<double> cosk1=(1. - sink1)/cosk0;
    //a0 = sqrt(b0+b0);
    //sink1= sqrt(sink1*(2.0-sink1))/sink0;
    sink0=(sink0-sink1)/(sink0+sink1);
    return(sink0.real()*sink0.real()+sink0.imag()*sink0.imag());
}

inline double multilayer::rf(vector<xrdata>::iterator pref0)
// reflectivity for a given angle
{

    // sink.at(0)=pref0->sink0;
    //cosk.at(0)=pref0->cosk0;
    int i=nl-1;
    //complex<double> a0=pref0->nsin2k,xs;
    complex<double> a0=pref0->nsin2k,xs;
    complex<double> r0(0,0),sink1,sink0;
    //xs=(a0 - nk[i])*(1.+nk[i]); // cos\theta/(1-nk) =1 - ( a0 +a0 nk -nk - nk^2)
    //sink1= sqrt(xs*(2.-xs)); //sin\theta = sqrt(1-cos^2\theta)
    //sink1 -= sink1*nk[i]; // sin\theta *(1- nk)
    sink1 = sqrt(a0 - nk[i]);
    complex<double> aki(1,0);
    //ak.at(--i)=complex<double>(1,0);
    while (--i>=0) {
        //xs=(a0 - nk[i])*(1.+nk[i]); // cos\theta/(1-nk) =1 - ( a0 +a0 nk -nk - nk^2)
        //sink0= sqrt(xs*(2.-xs)); //sin\theta = sqrt(1-cos^2\theta)
        //sink0 -= sink0*nk[i]; // sin\theta *(1- nk)
        sink0 = sqrt(a0 - nk[i]);
        xs=(sink0-sink1)/(sink0+sink1);
        sink1=sink0;
        r0 *= aki;
        r0 = (r0+xs)/(1.+xs*r0);
        aki=exp(kdz0*sink0);
    }
    return( r0.real()*r0.real()+r0.imag()*r0.imag());
}



void multilayer::genomerf(GARealGenome * g, double qmin, double qmax, int nq)
// output for forward calculation
{
    //if(g0.length() < nl+nl -3) {
    //cout<<"Incorrect genome size "<<g0.size()<<endl;
    //}
    //cout<<dz0<<endl;
    mkdensity(g);
    unsigned int j,jj=nl-1;
    double x=-1.5*dz0;
    double rnorm=1./n_bulk1.real(),inorm=1./n_bulk1.imag();
    ofstream out(fnrho.c_str());
    for (j=0;j<=jj;j++) {
        double r0;
        int j0= (int) ( x/dz1+0.5)-dside;
        if (j0<1) r0=0.;
        else r0=j0<glength-1?g->gene(j0):1.0; // the last gene is for bulk density
        r0 = rho_gene(r0).real()*rnorm;
        out<<x<<' '<<nk[j].real()*rnorm<<' '<<r0<<' '<<nk[j].imag()*inorm<<endl;
        x+=dz0;
    }
    out.close();
    out.open(fnrf.c_str());
    if (nq<1) nq=1;
    double q;
    if (qmax<qmin) {
        q=qmin;
        qmin=qmax;
        qmax=q;
    }
    q=qmin;
    double dq=(qmax-qmin)/nq;
    if (qmax==qmin) nq=0;
    //cout<<"qmin="<<qmin<<" qmax="<<qmax<<" nq="<<nq<<endl;
    vector<xrdata> ref2;
    int i=0;
    do {
        // cout<<b*180./M_PI<<" "<<a<<endl;
        // weighted using dyi[j], E(y) =  yi / dyi,
        double t=q*lambda/(4*M_PI);
        complex<double> a0=n_bulk0+t*t;
        double fresnel0=1./fresnel(t);
        ref2.push_back(xrdata(q,1.,1.,1.,1.,t,fresnel0,t,a0));
        q += dq;
    } while (i++<nq);


    vector<xrdata>::iterator pref2=ref2.begin();
    while (pref2 != ref2.end()) {
            double f=fresnel(pref2->sink0);
        out<<pref2->xi<<' '<<rf(pref2)/f<<' '<<f<<endl;
        pref2++;
    }
    out.close();
}
void multilayer::genomerf(GARealGenome * g)
// output for plotting
{
    //if(g0.length() < nl+nl -3) {
    //cout<<"Incorrect genome size "<<g0.size()<<endl;
    //}
    //cout<<dz0<<endl;
    mkdensity(g);
    unsigned int j,jj=nl-1;
    double x=-1.5*dz0;
    double rnorm=1./n_bulk1.real(),inorm=1./n_bulk1.imag();
    ofstream out(fnrho.c_str());
    for (j=0;j<=jj;j++) {
        double r0;
        int j0= (int) ( x/dz1+0.5)-dside;
        if (j0<1) r0=0.;
        else r0=j0<glength-1?g->gene(j0):1.0;
        r0 = rho_gene(r0).real()*rnorm;
        out<<x<<' '<<nk[j].real()*rnorm<<' '<<r0<<' '<<nk[j].imag()*inorm<<endl;
        x+=dz0;
    }
    out.close();
    //cout<<"generating "<< fnrf<<endl;
    out.open(fnrf.c_str());
    vector<xrdata>::iterator pref1=ref1.begin();
    while (pref1 != ref1.end()) {
     //   cout<<pref1->xi<<' '<<rf(pref1)<<' '<<fresnel(pref1->sink0)<<endl;
        out<<pref1->xi<<' '<<rf(pref1)/fresnel(pref1->sink0)<<endl;
        pref1++;
    }
    out.close();
}
complex<double> multilayer::rho_gene(double x)
//get density from gene
{
    if (x>1.) {
        double x0 =x-1;
        return( n_bulk1+(n_bulk12-n_bulk1)*x0);
    } else {
        return( n_bulk0+(n_bulk1-n_bulk0)*x);
    }
}

void multilayer::mkdensity(GARealGenome * g)
// generate density profile from genome
{
    double sigma0=exp(2.*g->gene(0)-1.0); //sigma
    double isigmas2=sqrt(0.5)/sigma0*dz1;
    double zrange2=6.*sigma0;
    double zrange=3.*sigma0;
    vector<complex<double> > nk0,dnk0;
    int i;
    zrange/=dz1;
    dside=(int)(zrange+0.5);
    nl=(int) ((slab+zrange2)/dz0)+2;
    zrange2/=dz1;
    nk.resize(nl);
    n_bulk1 =n_bulk10 * exp(0.5* (g->gene(g->size()-1) -1.)); // adjust bulk density
    nk.push_back(n_bulk0);
    nk0.push_back(n_bulk0);
    for (i=1;i<g->size()-1;i++) nk0.push_back(rho_gene(g->gene(i)));
    nk0.push_back(n_bulk1);
    for (i=1;i<nk0.size();i++) dnk0.push_back( 0.5*(nk0.at(i) - nk0.at(i-1)));
    complex<double> d0;
    int j;
    double z,idz1=dz0/dz1;
//#pragma omp parallel for private( i,j,z,zrange,isigmas2,d0 )
    for (i=0;i<nl;i++) {
        z=(i-0.5)*idz1;
        int il= (int) (z-zrange2)-1;
        if (il<0) il=0;
        int ir= (int) (z)+1;
        if (ir>dnk0.size()) ir=dnk0.size();
        d0=nk0.at(il);
        z -= il+zrange;
        for (j=il;j<ir;j++) {
            d0 += (1.+gsl_sf_erf (z*isigmas2))*dnk0.at(j);
            z -= 1.;
        }
        nk.at(i)=d0;
    }
    nk.at(nl-1)=n_bulk1;
    //debug
    /*
    ofstream out1("tmp.txt");
    for(i=0;i<nl;i++) out1<<i<<' '<<nk.at(i).real()/n_bulk1.real()<<endl;
    out1.close();
    */
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
    int nmax=ref0.size();
    //double d= g0.gene(0);
#pragma omp parallel for reduction(+:sy,sy2) private( a ) schedule(static)
    //while(pref0 !=ref0.end()){
    //for(vector<xrdata>::iterator pref0=ref0.begin(); pref0 !=ref0.end();pref0++){
    for (int i=0;i<nmax;i++) {
        //a = exp(d*x2i[j])* ( 1./6.*(rf(b - 0.005*M_PI/180.)+rf(b + 0.005*M_PI/180.))+4./6.*rf(b));
        //a = log(rf(pref0)/pref0->yi);
        a = log(rf(ref0.begin()+i)/ref0.at(i).yi);
        //cout<<thetai[j]<<" , "<<a<<endl;
        // weighted using dyi[j], E(y) =  yi / dyi,
        //               disabling weight
        //sy += a* pref0->wt;
        //sy2 += a*a* pref0->wt;
        sy += a* ref0.at(i).wt;
        sy2 += a*a* ref0.at(i).wt;
    }

    sy *= sdyi;
    g->y_factor=sy;
    sy2 *= sdyi;
    sy2 -= sy*sy;
    //a= sy*sy*(sy2 - sy*sy);
    //cout<<"ob: "<<a<<endl;
    return(sy2);
}

void multilayer::setbulk(double r0,double b0,double r1,double b1,double r12,double b12)
//set bulk phases times 2 to use for n^2
{
    n_bulk0= complex<double>(2.*2.8179e-5*r0*lambda*lambda/(2.*M_PI),-2.*b0);
    n_bulk10=n_bulk1=nk.at(nl-1)=complex<double>(2.*2.8179e-5*r1*lambda*lambda/(2.*M_PI),-2.*b1);
    n_bulk12=complex<double>(2.*2.8179e-5*r12*lambda*lambda/(2.*M_PI),-2.*b12);
}



#include <iostream>    // for cout, endl
#include <algorithm>   // for copy
#include <iterator>    // for ostream_iterator, xxx_inserter
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<cmath>
#include<sys/types.h>
#include<sys/wait.h>
#include<sys/stat.h>
#include<stdlib.h>
#include <grace_np.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;


#ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#  define EXIT_FAILURE -1
#endif

void my_error_function(const char *msg)
{
    fprintf(stderr, "library message: \"%s\"\n", msg);
}

class zprofile {
public:
    vector <double> xi, si, rri,rii;
    double dx,xmin,xmax,sigma,f0,f1,dx0;
    double zmin,zmax;
    double score;
    zprofile(){
    };
    void read( string s0);
};



void zprofile::read(string s0)
{
    ifstream in(s0.c_str());
    if( ! in.is_open()){
        cout<<"Can not open "<<s0<<endl;
        exit(0);
    }
    double lx,ly,lz,lw;
    in >> lx>>ly>>lz>>lw>>xmax;
    score=ly;
    xmin=zmin=lz;
    zmax=lw;
    while( in >> lx>>ly>>lz>>lw) {
        xi.push_back(lx);
        si.push_back(ly);
        rri.push_back(lz);
        rii.push_back(lw);
        if(xmin> lx - 3*lw) xmin=lx - 3.*lw;
        if(xmax< lx + 3*lw) xmax=lx + 3.*lw;

    }
    in.close();
    if ( xi.size() <2) {
        cout<<"Incorrect profile file "+s0<<endl;
        exit(0);
    }
    dx = (xmax - xmin)/512.;
    //dx= (xmax - xmin)/( xi.size() -1);
    //for(int i=0;i<xi.size();i++) cout<<xi.at(i)<<" "<<yi.at(i)<<endl;
}

int main(int argc, char *argv[])
{
    if(argc<4) return(0);

    GraceRegisterErrorFunction(my_error_function);
    string sfn(argv[1]+string("pop.dat"));
    string fnref(argv[1]);
    zprofile zf0;
    zf0.read(sfn);
    //for(int i=0;i<zf0.xi.size();i++) cout<<zf0.xi.at(i)<<" "<<zf0.yi.at(i)<<endl;
    string sfn2=fnref+".eps";
    string sid2(argv[2]);
    if (GraceOpenVA("gracebat", 16384, "-nosafe", "-noask",NULL)==-1){
    //if (GraceOpen(16384)==-1){
        cerr<<"Can't run Grace. \n";
        exit(EXIT_FAILURE);
    }
    GracePrintf("subtitle \" %s (\\f{Symbol}c\\f{}\\S2\\N= %g) \"",argv[3],zf0.score);
    cout<<"Objective= "<<zf0.score<<endl;
    double xmin,xmax;

    xmin= - 5.;
    xmax= zf0.xmax + 5.;
    double dx = (xmax - xmin)/512;
    cout<<xmin<<" "<<xmax<<endl;
    GracePrintf("page size 612, 792");
    GracePrintf("ARRANGE(2, 1, 0.15, 0.15,0.40)");
    //plotting rho
    GracePrintf("with g1");
    string fnrho=string(argv[1])+"rho.dat";
    ifstream inrho(fnrho.c_str());
    if( !inrho.is_open()){
        cout<<"Can not open "<<fnrho<<endl;
        exit(0);
    }
    for( unsigned int i1=0;i1<3;i1++){
        GracePrintf("s%d on",i1);
        GracePrintf("s%d symbol 0",i1);
        GracePrintf("s%d line linestyle 1",i1);
        GracePrintf("s%d line linewidth 2",i1);
    }
    double ymin,ymax,xi0,yi0,zi0;
    unsigned char xc;
    string linebuf;
    unsigned int il=0;
    vector <double> zi,rhoi;
    while (getline(inrho, linebuf)) {
        //cout<<il<<":";
        istringstream iss (linebuf);
        vector<double> va;
        std::copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(va));
        //            for(unsigned int jj=0;jj<va.size();jj++) cout<<' '<<va[jj];
        //           cout<<endl;
        if (va.size()>=2) {
            GracePrintf("s0 point %g,%g",va.at(0),va.at(1));
            if(va.size()>3) {
                GracePrintf("s1 point %g,%g",va.at(0),va.at(2));
                GracePrintf("s2 point %g,%g",va.at(0),va.at(3));
            }
            if(!il){
                xmax=xmin=va.at(0);
                ymax=ymin=va.at(1); //First line
            } else {
                if(xmin>va.at(0)) xmin=va.at(0);
                if(xmax<va.at(0)) xmax=va.at(0);
                if(ymin>va.at(1)) ymin=va.at(1);
                if(ymax<va.at(1)) ymax=va.at(1);
            }
        }
        zi.push_back(va.at(0));
        rhoi.push_back(va.at(1));
        il++;
    }
    inrho.close();
    /*
    {
            int ii1=0,ii2,ii3;
            double h0,h1;
    for(unsigned int ii=0;ii<zi.size();ii++){
    cout<<ii<<endl;
            if (ii1){
                    if(rhoi.at(ii)<0.5*ymax) {
                            ii3=ii;
                            h1= zi[ii-1] + (zi.at(ii) - zi[ii-1]) /(rhoi.at(ii)- rhoi[ii-1])*(0.5*ymax-rhoi[ii-1]);
                            cout<<"ii1 "<<ii<<endl;
                            break;
                    }
            } else {
                    if(rhoi.at(ii)>0.5*ymax) {
                            ii2=ii;
                            h0= zi[ii-1] + (zi.at(ii) - zi[ii-1]) /(rhoi.at(ii)- rhoi[ii-1])*(0.5*ymax-rhoi[ii-1]);
                            ii1=1;
                            cout<<"ii0 "<<ii<<endl;
                    }
            }
    }
    cout<<"FWHM = "<<h1-h0<<" angstrom\n";
    }
    */

    cout<<xmax<<" "<<ymax<<endl;
    if (xmax == xmin || ymax == ymin){
        cout<<"Incorrect format in "<<fnrho<<endl;
        exit(0);
    }


    GracePrintf("xaxis tick major %g",fabs(xmax - xmin) >40.?25.:5.);
    GracePrintf("world xmax %g",xmax);
    GracePrintf("world xmin %g",xmin);
    GracePrintf("world ymax %g",ymax*1.5);
    //GracePrintf("world ymax %g",(double) 15.);
    //GracePrintf("world ymin %g",ymin);

    //GracePrintf("autoscale");
    //GracePrintf("yaxis tick major %g",ymax>10.?(double)2:(double)1.);
    GracePrintf("autoticks");
    //cout<<"xmax= "<<xmax<<endl;
    GracePrintf("yaxis ticklabel char size 1.5");
    GracePrintf("xaxis ticklabel char size 1.5");
    GracePrintf("xaxis label \"\\+\\+z (\\f{Times-Roman}%c\\f{})\"",(unsigned char) 197);
    GracePrintf("yaxis label \"\\+\\+\\f{Symbol} r\\f{}(z)/\\f{Symbol}r\\f{}(\\f{Symbol}%c\\f{})\"",(unsigned char) 0xa5);
    GracePrintf("s2 hidden true");
    //GracePrintf("s1 hidden true");
    //wait(NULL);
    GracePrintf("redraw");
    GraceFlush();
    GracePrintf("with g0");
    //Reflectivity
    //string fn1("rsync -zuvt abraxas:"+spath+"/"+sfn+" .");
    //cout<<fn1<<endl;
    //system(fn1.c_str());
    string fnref0 = fnref + "ref.dat";
    ifstream in1(fnref0.c_str());
    if(!in1.is_open()){
        cout << "Can not open " << fnref0 <<endl;
        exit(0);
    }
    sfn=fnref+"rf.dat";
    ifstream in2(sfn.c_str());
    if(!in2.is_open()){
        cout << "Can not open " << sfn << endl;
        exit(0);
    }
    vector <double> xe,ye,dye,xf,yf;
    while (getline(in1, linebuf)) { // read in experimental data
        //cout<<il<<":";
        istringstream iss (linebuf);
        vector<double> va;
        std::copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(va));
        //            for(unsigned int jj=0;jj<va.size();jj++) cout<<' '<<va[jj];
        //           cout<<endl;
        if (va.size()>=3) {
            if(va.at(1)>=1.e-14) {
                xe.push_back(va.at(0));
                ye.push_back(va.at(1));
                dye.push_back(fabs(va.at(2)));
            }
        }
    }
    {//do we need to sqrt(errorbar)?
        unsigned int jj=0;
        for(unsigned int ii=0;ii<ye.size();ii++){
            if(dye.at(ii)<=0.3*ye.at(ii)*ye.at(ii) && dye.at(ii) < 1e-4*ye.at(ii)) jj++;
        }
        if( jj> (unsigned int) ((double) 0.3*ye.size())) {
            //yes, we need to
            cout<<"sqrt(dy)\n";
            for(unsigned int ii=0;ii<ye.size();ii++){
                dye.at(ii)=sqrt(dye.at(ii));
            }
        }
    }
    in1.close();
    while (getline(in2, linebuf)) { // read in fitting
        //cout<<il<<":";
        istringstream iss (linebuf);
        vector<double> va;
        std::copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(va));
        if (va.size()>=2) {
            xf.push_back(va.at(0));
            yf.push_back(va.at(1));
        }
    }
    in2.close();
    //begin interpolation of the fitting curve, make the values at experimental positions
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, xf.size());
    gsl_spline_init (spline, &(xf.at(0)), &(yf.at(0)), xf.size());
    ymax=0.;
    xmax=0.;
    for(unsigned int i2=0;i2<xe.size();i2++){
    if(xe.at(i2)>0.05){
        double a= dye.at(i2)/ye.at(i2);
        a= 1./(a*a);
        ymax += gsl_spline_eval (spline, xe.at(i2), acc)/ye.at(i2)*a;
        xmax += a;
	}
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    if(xmax == 0.) {
        cout<<"Divid by zero error, Check the input experimental data "<<fnref<<" now!\n";
    }
    ymax /= xmax;
    cout<<"Using "<<xmax<<", Rescale experimental data y *= "<<ymax<<endl;
    /*
    if(argc>=5) {
    	ymax=strtod(argv[4],NULL);
    }else {
    	if( xmax == 0.) ymax=1.;
    }
    */
    cout<<"Rescale experimental data y *= "<<ymax<<endl;
    for(unsigned int i2=0;i2<xe.size();i2++){
        ye.at(i2) *= ymax;
        dye.at(i2) *= ymax;
    }
    GracePrintf("s0 on");
    GracePrintf("s0 type xydy");
    GracePrintf("s0 errorbar size 0.15");
    GracePrintf("s1 on");
    GracePrintf("s1 type xy");
    ymax=ymin=ye.at(0);
    xmin=xe.at(0);
    xmax=xe.at(xe.size()-1);
    for(unsigned int i2=0;i2<xe.size();i2++){
        GracePrintf("g0.s0 point %g,%g",xe.at(i2),ye.at(i2));
        if(ymin> ye.at(i2)) ymin=ye.at(i2);
        if(ymax< ye.at(i2)) ymax=ye.at(i2);
        GracePrintf("g0.s0.y1[g0.s0.length-1] = %g",dye.at(i2));
    }
    for(unsigned int i2=0;i2<xf.size();i2++){
        GracePrintf("g0.s1 point %g,%g",xf.at(i2),yf.at(i2));
    }




    //GracePrintf("read xydy \"%s\"",fnref.c_str());
    //GracePrintf("s0.y1=sqrt(abs(s0.y))");
    //GracePrintf("read xy \"%s\"",sfn.c_str());
    //GracePrintf("yaxis  label \"\\+\\+|R\\s\\N\\+\\+(\\f{Symbol}q\\f{})|\"");
    GracePrintf("xaxis  ticklabel char size 1.5");
    GracePrintf("yaxis  ticklabel char size 1.5");
    GracePrintf("YAXES SCALE LOGARITHMIC");
    GracePrintf("autoticks");
    GracePrintf("yaxis  tick major 10");
    GracePrintf("yaxis  tick minor 9");
    GracePrintf("yaxis  ticklabel format power");
    GracePrintf("yaxis  ticklabel prec 0");
    GracePrintf("s0 line linestyle 0");
    GracePrintf("s1 line linestyle 1");
    GracePrintf("s1 line linewidth 2");
    GracePrintf("s0 symbol 1");
    GracePrintf("s0 symbol size 0.3");
    GracePrintf("world ymax %g",ymax*2.);
    //ymin=(double) 1.e-10;
    cout<<"xmax= "<<xmax<<"\tymin= "<<ymin<<endl;
    //GracePrintf("world ymin %g",0.25*ymin);
    GracePrintf("world ymin %g",ymin*0.5);
    GracePrintf("world xmax %g",0.2* ( (int) ( xmax/0.2) +1));
    GracePrintf("world xmin 0.");
    GracePrintf("xaxis tick major %g",xmax>1?0.25:0.1);
    //stringstream tsout;
    //for(int i=190;i<=200;i++){
    //tsout<<i<<'('<<(unsigned char) i<<')';
    //}
    //cout<<tsout.str()<<endl;
    GracePrintf("xaxis label \"\\+\\+q (%c\\S-1\\N\\+\\+)\"",(unsigned char) 197);
    GracePrintf("yaxis label \"\\+\\+R(q\\sz\\N\\+\\+)/R\\sF\\N\"");
    /*
    GracePrintf("yaxis  tick spec type both");
    GracePrintf("yaxis  tick spec 120");
    int i=0,j=1,k=0;
    double yc=j*pow(10.,i);
    do{
            if(j==1) {
                    GracePrintf("yaxis tick major %d,%g",k,yc);
                    if(i<0) {
                            if( (i/2)*2 == i) GracePrintf("yaxis ticklabel %d,\"10\\S%d\"",k,i);
                    }
                    else
                            GracePrintf("yaxis ticklabel %d,\"1\"",k);
                    k++;
                    i--;
                    j=9;

            } else {
                    GracePrintf("yaxis tick minor %d,%g",k++,yc);
                    j--;
            }
            yc=j*pow(10.,i);
    } while (yc>=ymin);
    */

    GracePrintf("device \"EPS\" dpi 300");
    GracePrintf("device \"EPS\" OP \"color,level2,bbox:tight\"");
    GracePrintf("hardcopy device \"EPS\"");
    GracePrintf("print to \"%s\"",sfn2.c_str());
    printf("print to \"%s\"\n",sfn2.c_str());
    GracePrintf("redraw");
    GracePrintf("print");
    //GracePrintf("exit");
    GraceFlush();

    if(GraceIsOpen()) GraceClose();
    waitpid(-1,NULL,0);
    /*
    for(int ii=0;ii<20;ii++){
    in1.open(sfn.c_str());
    if(in1.is_open()){
    in1.close();
    break;
    }
    sleep(1);
    }
    */
    //chmod the output file
    /*
    if(chmod(sfn2.c_str(),S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH )) {
    	cout<<"chmod "<<sfn2<<" error!\n";
    }
    */
    cout<<"Generating "<<sid2<<endl;
    pid_t pid = fork ();
    if (pid == 0)
    {//child
        execl("/usr/bin/convert","convert", sfn2.c_str(),sid2.c_str(),NULL);
    } else {
        if (pid<0) cout<<"fork() error!\n"; else waitpid(pid,NULL, 0);
    }
    if(chmod(sid2.c_str(),S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH )) {
        cout<<"chmod "<<sid2<<" error!\n";
    }
    return(0);
}

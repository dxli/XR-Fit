#pragma once
class xrdata{
public:
	double xi,yi,yin,dyn,dyi,rf,thetai,fresneli,wt,sink0;
	complex<double> nsin2k;
	xrdata(double xi0,double yin0,double dyn0,double yi0,double wt0,double thetai0,double fresneli0,double sink00,complex<double> nsin2k00){
		xi=xi0;
		yin=yin0;
		dyn=dyn0;
		yi=yi0;
		wt=wt0;
		thetai=thetai0;
		fresneli=fresneli0;
		sink0=sink00;
		nsin2k=nsin2k00;
	};
};


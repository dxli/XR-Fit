

#include "nph.h"
#include "nph_mysql.h"
#include "Task.h"
#include "elements.h"


// Multilayer reflection and transmission
// Parrett scheme, L. G. Parrett, Phys. Rev. 95(2), 359(1954)
// formulaa to use
//      n \cos(\theta) = \cos(\theta_0)
//
//      (E_{k,+} + E_{k,-} ) \cos(\theta_k) = ( E_{k+1,+}/a + E_{k+1,-} a) \cos( \theta_{k+1})
//      (E_{k,+} - E_{k,-} ) \sin(\theta_k) = ( E_{k+1,+}/a - E_{k+1,-} a) \sin( \theta_{k+1})
//      a = exp( i k_{k+1} d_{k+1})
//      k_{k+1} = sqrt( n^2 - \cos^2(\theta_0)) k_0

int main(int argc, char **argv)
{
    try
    {
        WriteLog("****************Starting NPH****************");
        run(atoi(argv[1]));
        WriteLog("********************DONE********************");
    }
    catch(exception &e)
    {
        string err = ("ERROR1: " + string(e.what()));
        WriteLog(err.c_str());
    }
    catch(int e)
    {
        ostringstream logStream;
        logStream << "ERROR2: " << e;
        WriteLog(logStream.str().c_str());
    }
}


void run(int taskId)
{

    ostringstream logStream;

    int i=total_slabs; // number of layers in Parrett reflectivity
    /*****************************************************************
     here is where to get the data input
     ******************************************************************/

    string data0;
    Task task0;
    getData(taskId, &data0, &task0);
    if(updateStarted(taskId, 1) == -1) {
        WriteLog("Error updating progress of task");
        exit(0);
    }
    double l0=E_to_l(task0.energy); //x-ray wavelength
    double dz0=task0.slab/(i-2);
    double dz1 =(( (unsigned int) 1)<<intp_factor)*dz0;
    multilayer ml0(i,l0,dz0,dz1);
    ostringstream stask; //making tmpfile names
    stask << taskId << "-";
    string fnref(("tmpout-" + stask.str()).c_str());
    ml0.fnpop=fnref+string("pop.dat");
    ml0.fnrf=fnref+string("rf.dat");
    string dl_folder("../htdocs/downloads/");
    stask.str("");
    stask<<dl_folder<<"ref-"<<taskId<<".txt";
    ml0.fnref=stask.str();
    stask.str("");
    ml0.fnrho=fnref+string("rho.dat");

    double rhoBulk0=task0.rho0;
    double rhoBulk1=task0.rho1;
    double rhoBulk12=task0.rho12;
    double betaBulk0=task0.beta0;
    double betaBulk1=task0.beta1;
    double betaBulk12=task0.beta12;
    //cout<<"rhoBulk0="<<rhoBulk0<<endl;
    //cout<<"rhoBulk1="<<rhoBulk1<<endl;
    unsigned int j;
    ml0.setbulk(rhoBulk0,betaBulk0,rhoBulk1,betaBulk1,rhoBulk12,betaBulk12);
    double roughness0=1.35+0.2*log(task0.slab/75.); //we use log(sigma0)
    unsigned int glength= (i>>intp_factor) ;
    ml0.glength=glength;

    ml0.qmin=task0.qmin;
    ml0.qmax=task0.qmax;
    ml0.readref(data0);

    if(ml0.ref0.size()<10) {
        WriteLog("Less than 10 data points read");
        updateDone(taskId, true);
        updateError(taskId, 1);
        exit(0);
    }
    data0.clear();



    // See if we've been given a seed to use (for testing purposes).  When you
    // specify a random seed, the evolution will be exactly the same each time
    // you use that seed number.
    char random_buf[257];


    string sr("/dev/urandom");
    struct timeval t_start;
    gettimeofday(&t_start,NULL);
    ifstream in1(sr.c_str());
    sr.resize(0);
    if(! in1.is_open())
    {
        WriteLog(("Can not open" + string(sr)).c_str());
        //cout<<"Can not open "<<sr<<endl;
        srandom(t_start.tv_usec);
    }
    else
    {
        in1.read(random_buf,256*sizeof(char));
        in1.close();
        initstate(t_start.tv_usec,random_buf,256);
        setstate(random_buf);
    }

    //cout<<"random: "<<random()<<endl;


    ofstream outfile;

    //  GARealAlleleSet alleles(0., MAX_VALUE);

    logStream.str("");
    //logStream << "gene: [ " << agene_min << ", " << agene_max << " ] ";
    WriteLog(logStream.str().c_str());
    logStream.str("");
    //cout<<"gene: [ "<<agene_min<<","<<agene_max<<"]\n";


    GARealGenome genome(1+glength , double(0.),double(2.)),
    genome1(1+glength , double(0.),double(2.)),
    genome2(1+glength , double(0.),double(2.));
    genome.randomize(MC_STEPSIZE);
    /*
        for(int jj=0;jj<genome.size();jj++) genome.gene(jj,(double) (g0[jj]* exp( -0.01*random()/(RAND_MAX+1.0))));
    */
    // dump the initial population to file

    //cout << "printing initial population to file..." << endl;

    genome1.copy(genome);
    int kk=0;
if(task0.progress>=1) {
    ifstream infile(ml0.fnpop.c_str());
    kk=1;

    if(infile.is_open())
    {
        double xi,yi,iyi;
        infile>>xi>>yi>>iyi;
        int jj=0;
        genome1.gene(jj++,iyi);

        while(jj< genome1.size())
        {
            if(infile>>xi>>yi>>iyi)
            {
                genome1.gene(jj,yi);
                jj++;
            }
            else
            {
                kk=0;
                break;
            }
        }
        infile.close();
    }
    else
        kk=0;
	}

    if(kk )
    {
        WriteLog(string("Using the saved genome in "+ml0.fnpop).c_str());
        //cout<<"Using the saved genome in "+ml0.fnpop<<endl;
    }
    else
    {
        WriteLog("Using a random genome");
        //cout<<"Using a random genome"<<endl;
        for(int jj=1;jj<genome1.size();jj++)
            genome1.gene(jj,(double) random()/(1.0+RAND_MAX)*2.0);
        genome1.gene(0,roughness0*0.2*(random()/(1.0+RAND_MAX)+4.5));
    }

    genome.copy(genome1);

    logStream << genome;
    WriteLog(logStream.str().c_str());
    logStream.str("");

    //cout<<genome<<endl;
    //ml0.genomerf(genome);
    //return 0 ;

    // cout<<"\n--magicalboundarystring\n";
    int ii=0;
    double score0,score1,score2,score20;
    //double score2_old;

    genome2.copy(genome);
    score20=score2=score0=ml0.objective(&genome);
    double mc_step=MC_STEPSIZE;
    unsigned acc0=0;
    unsigned int jmax=genome.size()<<2;

    /******************************************************************************

        Main loop!
        
    *******************************************************************************/
    gettimeofday(&t_start,NULL);
    unsigned int isteps=task0.progress;
    unsigned int isteps0=isteps;
    vector<double> dgene;
    dgene.resize(genome.size());
    for(int jj=0;jj<dgene.size();jj++) dgene.at(jj)=mc_step;
    while(isteps < task0.maxisteps)
    {
        for(unsigned int jj=0; jj<jmax;jj++)
        {
            int jg=(int) ( random()/(RAND_MAX+1.0)*(1.25-0.25*exp(-score0/0.05))*dgene.size());
	    if(jg>=dgene.size()) jg=0;
            double dx= dgene.at(jg)*random()/(RAND_MAX+1.0);
            double scorem,scorep;
            genome1.copy(genome);
            genome1.adjust(jg,dx);
            scorep=ml0.objective(&genome1)-score0;
            //genome1.copy(genome);
            genome1.adjust(jg,-dx - dx);
            scorem=ml0.objective(&genome1)-score0;
            if( scorem>0. && scorep>0. ){
                //double a0 = (scorem - scorep)/(2.*(scorem+scorep)); //second order
                genome1.adjust(jg,(3*scorem+scorep)/(2.*(scorem+scorep))*dx);
                if(score0<0.05) dgene.at(jg) *= scale_factord;
		else dgene.at(jg) *= 0.75;
            }else{
                if(scorem <0.) {
                    score0 += scorem;
                    genome.adjust(jg,-dx);
                    genome1.adjust(jg,-dx);
		    /*
                    genome.copy(genome1);
		    */
                }else{
                    score0 += scorep;
                    //genome.copy(genome1);
                    genome.adjust(jg,dx);
                    genome1.adjust(jg,3.*dx);
                }
                dgene.at(jg) *= scale_factori;
            }
            score1=ml0.objective(&genome1);
            if( score1< score0)
            {
                acc0++;
                //genome.copy(genome1);
                genome.gene(jg,genome1.gene(jg));
		genome.y_factor=genome1.y_factor;
                score0=score1;
            }
        }
        //double acpt=(double)acc0/jmax;
        if(score0<score20)
        {
            genome2.copy(genome);
	    score2=score20;
            score20=score0;
        }

        // cout<<isteps<<" acpt="<<acpt<<' '<<score0<<'('<<score2<<","<<score20<<") "<<mc_step<<endl;
        //if(acpt>0.3) mc_step *=1.2;
        //if( ml0.smooth_factor <0.55) ml0.smooth_factor +=0.02;
        //if(!acc0 || random()/(RAND_MAX+1.0) < score20)
        if(!acc0 || (score0>0.02 && (1-score20/score2)<0.05))
        {
            mc_step=MC_STEPSIZE;
            // cout<<genome<<endl;
            genome.randomize(mc_step);
            for(int ii=1;ii<dgene.size();ii++) dgene.at(ii)=mc_step;
            score0=ml0.objective(&genome);
            //cout<<"after randomize():\n"<<genome<<endl;
        }
        acc0=0;

        outfile.open(ml0.fnpop.c_str());
        int jj=0;
        double x=0.;
        outfile<<ii<<' '<<score20<<' '<<genome2.gene(jj++)<<endl;
        do
        {
            outfile<<x<<' '<<genome2.gene(jj)<<" 0.\n";
            x += ml0.dz1;
            jj++;
        }
        while(jj<genome2.size());

        outfile.close();


        isteps++;

        //generate image
        int status;
        pid_t pid;
        waitpid(-1,&status, WNOHANG);
        ml0.genomerf(&genome2);
        pid = fork ();
        if (pid == 0)
        {//child to generate image
            ostringstream sfnpng;
            sfnpng<<dl_folder<<"image-"<<taskId<<'-'<<isteps<<".png";
            string fnpng=sfnpng.str();
            pid_t pid2;
            pid2 = fork();
            if( pid2 == 0){
                timeval t_new;
                gettimeofday(&t_new,NULL);
                ostringstream logs1,ystr0;
                double t0=((t_new.tv_sec - t_start.tv_sec)+1e-6*(t_new.tv_usec - t_start.tv_usec))/(isteps-isteps0);
                logs1.str("");
                logs1<<fnpng<<" ("<<t0<<" s/Step)";
                WriteLog(logs1.str().c_str());
                cout<<logs1.str()<<endl;
                logs1.str("");
                logs1<<isteps<<'/'<<task0.maxisteps<<" ("<<t0<<" s/Step),("<<score0<<"),";
		ystr0<<exp(genome2.y_factor);
                //WriteLog(ystr0.str().c_str());
                //WriteLog("Plotting: ./plotmc1 ");
                //cerr<<"./plotmc1 "<<fnref.c_str()<<' '<<sid.str().c_str()<<endl;
               // cout<<"./plotmc1"<<"plotmc1"<<fnref.c_str()<<fnpng.c_str()<<logs1.str().c_str()<<ystr0.str().c_str()<<NULL;
                execl("./plotmc1","plotmc1",fnref.c_str(),fnpng.c_str(),logs1.str().c_str(),ystr0.str().c_str(),NULL);
            } else {
                if(pid2<0) {
                    WriteLog("fork error!");
                    exit(0);
                }
                waitpid(pid2,&status,0);
                chmod(fnpng.c_str(),S_IRUSR | S_IWUSR |S_IRGRP| S_IROTH);
                /*
                                     logStream.str("");
                                     if(errno0) {
                                         logStream<<"Unable to chmod "<<fnpng<<" : returns "<<errno0;
                                         WriteLog(logStream.str().c_str());
                                     }
                */
                updateProgress(taskId, isteps,((int) (100.*exp(-1+2.*genome.gene(0))+0.5))/100.);
		sfnpng.str("");
		sfnpng<<'-'<<taskId<<'-'<<isteps<<".txt";
		int icopy=0;
		string fns0("rf"),fns1(".dat");
		//cout<<"cp "<<fnref+fns0+fns1<<' '<<dl_folder+fns0+sfnpng.str()<<endl;
		icopy += file_copy(fnref+fns0+fns1,dl_folder+fns0+sfnpng.str());
		fns0=string("rho");
		//cout<<"cp "<<fnref+fns0+fns1<<' '<<dl_folder+fns0+sfnpng.str()<<endl;
		icopy += file_copy(fnref+fns0+fns1,dl_folder+fns0+sfnpng.str());
		//cout<<"icopy: "<<icopy<<endl;
		if(icopy) {
			WriteLog("ERROR: data file copying");
		}
		}
            exit(0);
            } else
        {
            if (pid<0)
            {
                WriteLog("fork() error!");
            }
            else waitpid(-1,&status, WNOHANG);
        }


        //output image

        //check "started"
        if(getStarted(taskId)<0 ) {
            cout<<"Stop Task "<<taskId<<endl;
            break;
        }

    }
    updateDone(taskId, true);
}


//reads data in from database and returns the tasks id
void getData(int id, string *data0, Task *t)
{
    string d;

    if(id > 0)
    {
        getTask(t, id);

        istringstream str(t->data.c_str());


        int len = t->data.size();
        data0->resize(len);
        str.read(&(*data0)[0], len);
    }
    else
    {
        WriteLog("getData - Error Retrieving Data");
    }
}

void read_in_data(string * data0, char ** argv, int argc)
{
    string fname = argv[1];
    int len = atoi(argv[2]);
    data0->resize(len);
    ifstream dfile(fname.c_str());
    dfile.read(&(*data0)[0], len);
    dfile.close();
}

ostream& operator << (ostream& os, GARealGenome genome)
{
    // Output the atom position to the stream in the format n/d
    for(int jj=0;jj<genome.size();jj++) os<<' '<<genome.gene(jj);
    return os;
}

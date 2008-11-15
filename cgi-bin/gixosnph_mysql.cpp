#include "nph_mysql.h"
#include "elements.h"


int getNextTaskId()
//searching for a task to start
{
    int id = -1,progress;
    ostringstream strbuf;
    strbuf << "SELECT id,progress FROM tasks WHERE started=0 and running=0 ORDER BY date ASC LIMIT 1";
    
    Connection con(use_exceptions);
    
    try
    {
	Query query = con.query();
	con.connect(DATABASE, HOST, USER, PASSWORD);
	query << strbuf.str();
 StoreQueryResult res = query.store();
	cout<<strbuf.str()<<endl;
	
        if (res && res.num_rows() > 0)
	{
		mysqlpp::Row row;
		row = res.at(0);
		id = row["id"];
		progress = row["progress"];
	}

    }
    catch (const BadQuery& er)
    {    
    // Handle any query errors
        cerr << "getNextTaskId - Query error: " << er.what() << endl;
	    id = -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "getNextTaskId - Conversion error: " << er.what() << endl <<
			    "\tretrieved data size: " << er.retrieved <<
			    ", actual size: " << er.actual_size << endl;
	    id = -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "getNextTaskId - Error: " << er.what() << endl;
	    id = -1;
    }
    if(id>=1 && ! (progress >=1)) { //first time run, convert formula/density
    try{
    	strbuf.str("");
    strbuf << "SELECT energy,formula0,formula1,formula12,rhoin0,rhoin1,rhoin12 FROM tasks WHERE id="<<id;
    	Query query = con.query();
	con.connect(DATABASE, HOST, USER, PASSWORD);
	query << strbuf.str();
 StoreQueryResult res = query.store();
	cout<<strbuf.str()<<endl;
	mysqlpp::Row row;
		row = res.at(0);
		string formula0(row["formula0"]), formula1(row["formula1"]), formula12(row["formula12"]);
		double energy=row["energy"];
		double rhoin0=row["rhoin0"];
		double rhoin1=row["rhoin1"];
		double rhoin12=row["rhoin12"];
		double rho0=0.,rho1=0.,rho12=0.;
		double beta0=0.,beta1=0.,beta12=0.;
		cout<<formula0.size()<<endl;
		cout<<"formula0= "<<formula0<<endl;
		cout<<"formula1= "<<formula1<<endl;
		cout<<"formula12= "<<formula12<<endl;
		cout<<"rhoin0= "<<rhoin0<<endl;
		cout<<"rhoin1= "<<rhoin1<<endl;
		cout<<"rhoin12= "<<rhoin12<<endl;
		compound::compound cmpd0(formula0,energy,rhoin0), cmpd1(formula1,energy,rhoin1), cmpd12(formula12,energy,rhoin12);
		strbuf.str("");
        strbuf << "UPDATE tasks SET rho0="<<cmpd0.rho_el<<",rho1="<<cmpd1.rho_el<<",rho12="<<cmpd12.rho_el<<",beta0="<<cmpd0.beta<<",beta1="<<cmpd1.beta<<",beta12="<<cmpd12.beta<<" WHERE id="<< id;
        //cout<<strbuf.str()<<endl;
        query.exec(strbuf.str());
	cout<<"rho0="<<cmpd0.rho_el<<"\tbeta="<<cmpd0.beta<<endl;
	cout<<"rho1="<<cmpd1.rho_el<<"\tbeta="<<cmpd1.beta<<endl;
	cout<<"rho12(Maximum)="<<cmpd12.rho_el<<"\tbeta="<<cmpd12.beta<<endl;
    }
    
       catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "updateDone - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "updateDone - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "updateDone - Error: " << er.what() << endl;
        return -1;
    }
}
 
    return id;
}


//allocate memory for Task *t before calling the function
void getTask(Task * t, int id)
{
    ostringstream strbuf;
    strbuf << "SELECT * FROM tasks WHERE id=" << id;
    cout<<strbuf.str()<<endl;
    
    Connection con(use_exceptions);
    
    try
    {
        Query query = con.query();
        con.connect(DATABASE, HOST, USER, PASSWORD);
        query << strbuf.str();
        StoreQueryResult res = query.store();
    
        if (res)
        {
            mysqlpp::Row row;
            row = res.at(0);
            t->id = row["id"];
cout<<row["id"]<<endl;
            t->progress = row["progress"];
            t->maxisteps = row["maxisteps"];
            t->date = string(row["date"]);
            t->email = string(row["email"]);
            t->data = string(row["data"]);
	    /*
            t->formula0 = string(row["formula0"]);
            t->formula1 = string(row["formula1"]);
            t->formula12 = string(row["formula12"]);
	    */
            t->started = row["started"];
            t->running = row["running"];
            t->done = row["done"];
            t->error = row["error"];
            t->rho0 = row["rho0"];
            t->rho1 = row["rho1"];
            t->rho12 = row["rho12"];
            t->beta0 = row["beta0"];
            t->beta1 = row["beta1"];
            t->beta12 = row["beta12"];
            t->slab = row["slab"];
            t->energy = row["energy"];
            t->qmin = row["qmin"];
            t->qmax = row["qmax"];
            t->nslab = row["nslab"];
            t->dth = row["dth"];
            t->alpha0 = row["alpha0"];
            t->gamma = row["gamma"];
            t->sigmad = row["sigmad"];
	    /*
            t->rhoin0 = row["rhoin0"];
            t->rhoin1 = row["rhoin1"];
            t->rhoin12 = row["rhoin12"];
	    */
        }
    }
    catch (const BadQuery& er)
    {    
    // Handle any query errors
        cerr << "getTask - Query error: " << er.what() << endl;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "getTask - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "getTask - Error: " << er.what() << endl;
    }
}



int updateProgress(int id, int progress,double sigma0)
{
    Connection con(use_exceptions);
    try
    {
        ostringstream strbuf;
        unsigned int i = 0;
        con.connect(DATABASE, HOST, USER, PASSWORD);
        Query query = con.query();
        strbuf << "UPDATE tasks SET progress=" << progress << ",sigma0="<<sigma0<<" WHERE id=" << id;
        query.exec(strbuf.str());
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "updateProgress - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "updateProgress - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "updateProgress - Error: " << er.what() << endl;
        return -1;
    }

    return 0;
}

int updateStarted(int id, int started)
{
    Connection con(use_exceptions);
    try
    {
        ostringstream strbuf;
        unsigned int i = 0;
        con.connect(DATABASE, HOST, USER, PASSWORD);
        Query query = con.query();
        strbuf << "UPDATE tasks SET started=" << started << " WHERE id=" << id;
        query.exec(strbuf.str());
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "updateStarted - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "updateStarted - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "updateStarted - Error: " << er.what() << endl;
        return -1;
    }

    return 0;
}
     
     
int updateDone(int id, bool done)
{
    Connection con(use_exceptions);
    try
    {
        ostringstream strbuf;
        unsigned int i = 0;
        con.connect(DATABASE, HOST, USER, PASSWORD);
        Query query = con.query();
        strbuf << "UPDATE tasks SET done=" << done << ",started=-1,running=0 WHERE id=" << id;
        query.exec(strbuf.str());
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "updateDone - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "updateDone - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "updateDone - Error: " << er.what() << endl;
        return -1;
    }

    return 0;
}        


 
int updateRunning(int id, int running)
{
    Connection con(use_exceptions);
    try
    {
        ostringstream strbuf;
        unsigned int i = 0;
        con.connect(DATABASE, HOST, USER, PASSWORD);
        Query query = con.query();
        strbuf << "UPDATE tasks SET running="<<running<<" WHERE id=" << id;
        query.exec(strbuf.str());
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "updateDone - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "updateDone - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "updateDone - Error: " << er.what() << endl;
        return -1;
    }

    return 0;
}        

 
int updateError(int id, int error)
{
    Connection con(use_exceptions);
    try
    {
        ostringstream strbuf;
        unsigned int i = 0;
        con.connect(DATABASE, HOST, USER, PASSWORD);
        Query query = con.query();
        strbuf << "UPDATE tasks SET error="<<error<<" WHERE id=" << id;
        query.exec(strbuf.str());
	strbuf.str("");
        strbuf << "UPDATE tasks SET done=0 WHERE id=" << id;
        query.exec(strbuf.str());
	strbuf.str("");
        strbuf << "UPDATE tasks SET started=-1 WHERE id=" << id;
        query.exec(strbuf.str());
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "updateDone - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "updateDone - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "updateDone - Error: " << er.what() << endl;
        return -1;
    }

    return 0;
}        
 
int getStarted(int id)
{
	int started;
    Connection con(use_exceptions);
    try{
    	ostringstream strbuf;
 	Query query = con.query();
        con.connect(DATABASE, HOST, USER, PASSWORD);
        strbuf << "select started from tasks WHERE id=" << id;
        query << strbuf.str();
        StoreQueryResult res = query.store();
    
        if (res)
        {
            mysqlpp::Row row;
            row = res.at(0);
            started = row["started"];
	    }
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "getStarted - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "getStarted - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "getStarted - Error: " << er.what() << endl;
        return -1;
    }

    return started;
}        




 
int getRunning(int id)
{
	int running;
    Connection con(use_exceptions);
    try{
    	ostringstream strbuf;
 	Query query = con.query();
        con.connect(DATABASE, HOST, USER, PASSWORD);
        strbuf << "select running from tasks WHERE id=" << id;
        query << strbuf.str();
        StoreQueryResult res = query.store();
    
        if (res)
        {
            mysqlpp::Row row;
            row = res.at(0);
            running = row["running"];
	    }
    }
    catch (const BadQuery& er)
    {
    // Handle any query errors
        cerr << "getRunning - Query error: " << er.what() << endl;
        return -1;
    }
    catch (const BadConversion& er)
    {
    // Handle bad conversions
        cerr << "getStarted - Conversion error: " << er.what() << endl <<
                "\tretrieved data size: " << er.retrieved <<
                ", actual size: " << er.actual_size << endl;
        return -1;
    }
    catch (const Exception& er)
    {
    // Catch-all for any other MySQL++ exceptions
        cerr << "getRunning - Error: " << er.what() << endl;
        return -1;
    }

    return running;
}        





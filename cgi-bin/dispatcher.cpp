#include "Task.h"
#include "nph_mysql.h"
#include "logUtility.h"


int main(int argc, char **argv)
{


    int id = -1,id0=-1;
    pid_t pid;

    while(true)
    {
        id = getNextTaskId();
        if(id>0 ) cout<<"next id="<<id<<endl;
        id0=id;
        if(id > 0)
        {
            pid = fork();
            if(pid == 0)
            {
                cout << "Starting task " << id << endl << endl;
                WriteLog("-----------------Starting------------------");
                ostringstream sid;
                sid << id;
                execl("./nph", "nph", sid.str().c_str(), NULL);
            }
            else
            {
                ostringstream log;
                log << "Started nph with pid " << pid;
                WriteLog(log.str().c_str());
                updateRunning(id,1);
		       int status;
        waitpid(pid,&status,0);
        int re0=getRunning(id);
	if(re0) {
		log.str("");
                log <<"Task "<<id<<" exited abnormally. Setting running=0,started=0 for auto launching";
                WriteLog(log.str().c_str());
		updateRunning(id,0);
		updateStarted(id,0);
	}
                id0=-1;
        }

        }
        else
        {
            cout << "no new tasks found. sleep(10)" << endl << endl;
	    sleep(10);
        }
        //cout << "sleeping 10 seconds" << endl << endl;
        //sleep(10);
     }
    return 0;
}

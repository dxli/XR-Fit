#include <pthread.h>


#include "nph.h"
#include "logUtility.h"
#include <stdio.h>
#include <exception>
#include <string>
#include <iostream>

using namespace std;


int read_form(string *d0);
void form_error(char *s0);


int main(int argc, char *argv[])
{
    
    try
    {
        WriteLog("-----------------Starting------------------");

        WriteLog("-----------------Read Form-----------------");

        string data;
        read_form(&data);

        int len = data.size();

        int pid = fork();

        if(pid != 0)
        {
            srand(time(0));
            int a = rand();
            ostringstream strStream;
            strStream << a << ".txt";
            ofstream dfile(strStream.str().c_str());
            dfile.write(data.c_str(), len);
            dfile.close();

            ostringstream slen;

            slen << len;

            WriteLog("-------------Staring Executable------------");

            execl("./nph", "nph", strStream.str().c_str(), slen.str().c_str(), NULL);
        }
        else
        {
            WriteLog("----------------Leaving CGI----------------");
        }
    
    }
    catch(exception &e)
    {
        string err = ("ERROR: " + string(e.what()));
        WriteLog(err.c_str());
    }
    catch(int e)
    {
        ostringstream logStream;
        logStream << "ERROR: " << e;
        WriteLog(logStream.str().c_str());
    }
    return 0;
}


int read_form(string *d0)
{
    char *lenstr;
    string data0;
    lenstr = getenv("CONTENT_LENGTH");

    if(lenstr==NULL)
    {
        WriteLog("CONTENT_LENGTH: 0 -- exiting");
        exit(0);
    }
    
    WriteLog( (string("CONTENT_LENGTH: ") + string(lenstr)).c_str() );
    
    
    
    int data_len;
    if( sscanf(lenstr,"%d",&data_len)!=1)
        exit(0);
    
    if( data_len<=0 || data_len>=20000)
        exit(0);
    
    data0.resize(data_len);
    
    cin.read(&(data0[0]),data_len);
    
    string s0("filename=\"");
    
    int i=data0.find(s0,0)+s0.size();
    
    if(i==string::npos)
        WriteLog("Invalid input");
    
    if(data0.at(i) == '"' )
    {
        // no file search for poste
        s0=string("name=\"pasted\"");
        i=data0.find(s0,0)+s0.size()+2;
    }
    else
    {
        //read uploaded
        i=data0.find("\n",i)+1;
        i=data0.find("\n",i)+1;
    }
    
    int j=data0.find("-----------------------------",i);
    
    *d0=data0.substr(i,j-i);
}





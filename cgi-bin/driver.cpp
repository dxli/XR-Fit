#include <pthread.h>


#include "nph.h"
#include <stdio.h>
#include <exception>
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    
    try
    {
        WriteLog("****************Starting NPH****************");
        driver(argc, argv);
        WriteLog("********************DONE********************");
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


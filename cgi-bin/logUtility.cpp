#include "logUtility.h"

void WriteLog(const char * s)
{
	ofstream log("nph.log", ios::app);
	time_t rawtime;
	time ( &rawtime );
	string t = ctime(&rawtime);
	t[t.length()-1] = '\0';
	log << "[ " <<  t << " ] : " << s << endl;
	log.close();
}




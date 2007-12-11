#include "nph.h"
using namespace std;
int file_copy(string fn_in,string fn_out)
//copy a file
{
    ifstream in0(fn_in.c_str());
    if(!in0.is_open()) return(1);
    //cout<<fn_in<<" opened\n";
    filebuf *pbuf=in0.rdbuf();
    size_t size=pbuf->pubseekoff (0,ios::end,ios::in);
    pbuf->pubseekpos (0,ios::in);
    if(size>=0x20000) {
        in0.close();
        return(3);
    }
    string buffer;
    //cout<<"file size="<<size<<endl;
    buffer.resize(size);
    pbuf->sgetn (&(buffer[0]),size);
    in0.close();
    ofstream out0(fn_out.c_str());
    if(!out0.is_open()) return(7);
    //cout<<fn_out<<" opened\n";
    out0.write(&(buffer[0]),size);
    out0.close();
    chmod(fn_out.c_str(),S_IRUSR | S_IWUSR |S_IRGRP| S_IROTH);
    return 0;
}

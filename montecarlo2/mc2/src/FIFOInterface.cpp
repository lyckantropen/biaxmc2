#include "FIFOInterface.h"

int FIFOInterface::ref = 0 ;


int FIFOInterface::Read(std::string &s) {
    char buff[255];
    int l = read(in,(void*)buff,255);
    if(l>0) s = std::string(buff);
    return l;
}

void FIFOInterface::Write(const char *cmd){
    Write(std::string(cmd));
}

void FIFOInterface::Write(const int &i){
    std::stringstream s;
    s << i << std::endl;
    Write(s.str());
}

void FIFOInterface::Write(const std::string &cmd){
    int out = open(out_name.string().c_str(),O_WRONLY|O_NONBLOCK);
    if(write(out,(const void*)(cmd.c_str()),cmd.length())==-1)
        std::cout << std::strerror(errno) << std::endl;
    fsync(out);
    close(out);
}

FIFOInterface::FIFOInterface(const std::string &bbase)
{
    in_name=fs::path(bbase)/"in";
    out_name=fs::path(bbase)/"out";
    if(!fs::exists(bbase))
        fs::create_directories(bbase);


    mkfifo(in_name.string().c_str(),0666);
    mkfifo(out_name.string().c_str(),0666);
    in = open(in_name.string().c_str(),O_RDONLY|O_NONBLOCK);

}

FIFOInterface::FIFOInterface(const FIFOInterface &s){
    in = s.in;
    in_name = s.in_name;
    out_name = s.out_name;
    ref++;
}

FIFOInterface::~FIFOInterface(){
    ref--;
    if(ref==0){
        close(in);
        fs::remove(in_name);
        fs::remove(out_name);
    }
}

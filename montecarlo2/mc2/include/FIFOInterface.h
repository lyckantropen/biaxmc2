/* 
 * File:   FIFOInterface.h
 * Author: karol
 *
 * Created on 9 wrzesień 2010, 12:41
 */

#ifndef _FIFOINTERFACE_H
#define	_FIFOINTERFACE_H

#include "stringsaveload.h"
#include "boost.h"
#include "std.h"
#include "sys/types.h"
#include "sys/stat.h"
#include "fcntl.h"
#include "unistd.h"
#include <cerrno>
#include <cstdio>
#include <cstring>

class FIFOInterface {
    //std::ifstream in;
    //std::ofstream out;

    static int ref;
    int in;
    fs::path in_name;
    fs::path out_name;
protected:
    int Read(std::string & s) {
        char buff[255];
        int l = read(in,(void*)buff,255);
        if(l>0) s = std::string(buff);
        return l;
    }
    template<class item_t>
    void Write(item_t & item){
        int out = open(out_name.string().c_str(),O_WRONLY|O_NONBLOCK);
        std::string s = boostbase::serialsave<item_t>(item);
        int w = write(out,(const void*)(s.c_str()),s.length());
        fsync(out);
        if(w<s.length())
            std::cout << w << "<" << s.length() << std::endl;
        close(out);
    }
    void Write(const char * cmd){
        Write(std::string(cmd));
    }
    void Write(const int & i){
        std::stringstream s;
        s << i << std::endl;
        Write(s.str());
    }
    void Write(const std::string & cmd){
        int out = open(out_name.string().c_str(),O_WRONLY|O_NONBLOCK);
        if(write(out,(const void*)(cmd.c_str()),cmd.length())==-1)
            std::cout << std::strerror(errno) << std::endl;
        fsync(out);
        close(out);
    }
    template<class item_t>
    item_t Restore(const std::string & s){
        return boostbase::serialload<item_t>(s);
    }
public:
    FIFOInterface(const fs::path & base):
    in_name(base/"in"),
    out_name(base/"out")
    {
        if(!fs::exists(base))
            fs::create_directories(base);
        mkfifo(in_name.string().c_str(),0666);
        mkfifo(out_name.string().c_str(),0666);
        in = open(in_name.string().c_str(),O_RDONLY|O_NONBLOCK);
        
    }
    FIFOInterface(const FIFOInterface & s){
        in = s.in;
        in_name = s.in_name;
        out_name = s.out_name;
        ref++;
    }
    
    ~FIFOInterface(){
        ref--;
        if(ref==0){
            close(in);
            fs::remove(in_name);
            fs::remove(out_name);
        }
    }
     

};



#endif	/* _FIFOINTERFACE_H */


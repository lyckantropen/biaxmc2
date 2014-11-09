/*
 * File:   FIFOInterface.h
 * Author: karol
 *
 * Created on 9 wrzesień 2010, 12:41
 */

#ifndef _FIFOINTERFACE_H
#define _FIFOINTERFACE_H

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

/// CURRENTLY UNUSED interface prototype
/// opens two FIFO file objects which provide a way to communicate
/// with the application through commandline or other applications
class FIFOInterface
{
    //std::ifstream in;
    //std::ofstream out;

    static int ref;
    int in;
    fs::path in_name;
    fs::path out_name;
protected:
    int Read(std::string & s);
    template<class item_t>
    void Write(item_t & item)
    {
        int out = open(out_name.string().c_str(), O_WRONLY | O_NONBLOCK);
        std::string s = boostbase::serialsave<item_t>(item);
        int w = write(out, (const void*)(s.c_str()), s.length());
        fsync(out);
        if(w < s.length())
            std::cout << w << "<" << s.length() << std::endl;
        close(out);
    }
    void Write(const char * cmd);
    void Write(const int & i);
    void Write(const std::string & cmd);
    template<class item_t>
    item_t Restore(const std::string & s)
    {
        return boostbase::serialload<item_t>(s);
    }
public:
    FIFOInterface(const std::string & bbase);
    FIFOInterface(const FIFOInterface & s);

    ~FIFOInterface();


};



#endif  /* _FIFOINTERFACE_H */


/* 
 * File:   ILoggable.h
 * Author: karol
 *
 * Created on 3 grudzień 2009, 12:30
 */

#ifndef _ILOGGABLE_H
#define	_ILOGGABLE_H

#include "std.h"
#include "sys/types.h"
#include "unistd.h"
#include "omp.h"

class ILoggable {
protected:
    std::stringstream log;
    io::null_sink	null_sink;
    io::stream<io::null_sink>	null;
    std::ostream * stream;
    ILoggable():null(null_sink)
	{
        //stream = &log;
	stream = &null;
    }

public:

    ILoggable(const ILoggable & s){
        log << s.log.str();
        if(s.stream==&s.log)
            stream=&log;
        else
            stream=const_cast<std::ostream*>(s.stream);
    }
    const ILoggable & operator=(const ILoggable & s){
        log << s.log.str();
        if(s.stream==&s.log)
            stream=&log;
        else
            stream=const_cast<std::ostream*>(s.stream);
        return *this;
    }
    std::string GetLog() const {
        return log.str();
    }
    std::string Print() {
        return log.str();
    }
    virtual std::ostream & Log() {
        return *stream;
    }
    virtual void SetStream(std::ostream * os){
        if(os!=NULL)
            stream = os;
        else
            stream = &log;d
    }
    virtual void SetFile(const std::string & f){
        std::stringstream s;
        s << getpid() << "_" << omp_get_thread_num() << "_" << f;
        stream = new std::ofstream(s.str().c_str());
    }
};

#endif	/* _ILOGGABLE_H */


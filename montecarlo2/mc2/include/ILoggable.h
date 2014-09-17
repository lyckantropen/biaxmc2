/*
 * File:   ILoggable.h
 * Author: karol
 *
 * Created on 3 grudzień 2009, 12:30
 */

#ifndef _ILOGGABLE_H
#define _ILOGGABLE_H

#include "std.h"
#include "sys/types.h"
#include "unistd.h"
#include "boost.h"

///
/// Creates a flexible logging mechanism for the
/// child class. The log can be either a file or a standard
/// output stream (such as std::cout).
///
///
template<class Holder>
class loggable_proxy
{
    std::mutex access;
    std::function<void(const std::string&)> f_print;
    bool pass;
public:
    loggable_proxy(Holder * _h, void (Holder::*f)(const std::string &))
    {
        pass = false;
        using std::placeholders::_1;
        f_print = std::bind(f,_h,_1);
    }
    template<class val_t>
    loggable_proxy<Holder> & operator<<(const val_t & s)
    {
        if(!pass)
        {
            std::stringstream str;
            str << s;
            f_print(str.str());
        }
        return * this;
    }
    loggable_proxy<Holder> & operator<<( std::ostream &(ff)(std::ostream &) )
    {
        if(!pass)
        {
            std::stringstream str;
            ff(str);
            f_print(str.str());
        }
        return * this;
    }
    void SetPass(bool p)
    {
        pass = p;
    }
};

class ILoggable
{
//    std::mutex print_mutex;

//    void print_thread(std::string s);

    loggable_proxy<ILoggable> proxy;

    std::stringstream log;
    io::null_sink   null_sink;
    io::stream<io::null_sink>   null;
    std::ostream * stream;
    bool internal;

public:

    ILoggable();
    ILoggable(const ILoggable & s);
    virtual ~ILoggable();
    const ILoggable & operator=(const ILoggable & s);
    std::string GetInternalLog() const;
    void Print(const std::string & s);

    /// Use this function as in Log() << "output";
    virtual loggable_proxy<ILoggable> &Log();
    void SetStream(std::ostream * os);
    void SetStream(ILoggable * other);
    void SetFile(const std::string & f);

    void SetInternal(bool i);
};

#endif  /* _ILOGGABLE_H */


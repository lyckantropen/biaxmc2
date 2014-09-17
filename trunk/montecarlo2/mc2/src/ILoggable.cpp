#include "ILoggable.h"
#ifndef OMP_DEBUG
#include <omp.h>
#else
#include "omp_debug.h"
#endif

//void ILoggable::print_thread(std::string s)
//{
//    std::lock_guard<std::mutex> lock(print_mutex);
//    (*stream) << s;
//}

ILoggable::ILoggable():
    proxy(this,&ILoggable::Print),
    null(null_sink),
    internal(false),
    stream(&null)
{
    proxy.SetPass(true);
}

ILoggable::ILoggable(const ILoggable &s):
    proxy(this,&ILoggable::Print),
    null(null_sink),
    internal(s.internal)
{
    if(s.stream==&s.null)
    {
        stream=&null;
        proxy.SetPass(true);
    }
    else
        stream=const_cast<std::ostream*>(s.stream);
}

ILoggable::~ILoggable()
{
}

const ILoggable &ILoggable::operator=(const ILoggable &s)
{
    log << s.log.str();
    internal = s.internal;
    if(s.stream==&s.null)
    {
        stream=&null;
        proxy.SetPass(true);
    }
    else
        stream=const_cast<std::ostream*>(s.stream);
    return *this;
}

std::string ILoggable::GetInternalLog() const
{
    return log.str();
}

void ILoggable::Print(const std::string &s)
{
//    std::thread t = std::thread::thread(std::bind(&ILoggable::print_thread,this,s));
//    t.detach();
    (*stream) << s;
    // TODO: duplicate log
    if(internal)
        log << s;
}

loggable_proxy<ILoggable> &ILoggable::Log() {
    return proxy;
}

void ILoggable::SetStream(std::ostream *os){
    if(os!=NULL)
    {
        stream = os;
        proxy.SetPass(false);
    }
    else
    {
        stream = &null;
        proxy.SetPass(true);
    }
}

void ILoggable::SetStream(ILoggable *other)
{
    SetStream(other->stream);
}

void ILoggable::SetFile(const std::string &f){
    std::stringstream s;
    s << getpid() << "_" << omp_get_thread_num() << "_" << f;
    stream = new std::ofstream(s.str().c_str());
}

void ILoggable::SetInternal(bool i)
{
    internal = i;
}

/* 
 * File:   std.h
 * Author: karol
 *
 * Created on 4 listopad 2009, 17:45
 */

#ifndef _STD_H
#define	_STD_H

#include <exception>
#include <iostream>
#include <sstream>
#include <iosfwd>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <utility>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <valarray>
#include <memory>
#include <type_traits>
#include <thread>
#include <functional>
#include <mutex>


using std::unique_ptr;
using std::weak_ptr;
using std::shared_ptr;
using std::make_shared;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;
using std::remove_const;

#if (__cplusplus <= 201103L)

/// a trivial definition of make_unique, since it is absent from C++11
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#endif

/// sleep m milliseconds
#ifndef SLEEP_MS
#define SLEEP_MS(m) (std::this_thread::sleep_for(std::chrono::milliseconds(m)))
#endif

#endif	/* _STD_H */


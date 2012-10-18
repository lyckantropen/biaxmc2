/* 
 * File:   serialhash.h
 * Author: karol
 *
 * Created on 8 listopad 2009, 14:08
 */

#ifndef _SERIALHASH_H
#define	_SERIALHASH_H

#include "boost.h"
#include "serializer.h"
#include "stringsaveload.h"
#include "std.h"
#include "hl_md5wrapper.h"

namespace boostbase {

    template<class item_t>
    std::string md5gen(item_t & item) {
        hashwrapper * md5 = new md5wrapper();
        std::string content = serialsave<item_t > (item);
        std::string hash = md5->getHashFromString(content);
        delete md5;
        return hash;
    }

    template<class item_t>
    bool md5check(item_t & item, const std::string & hash) {
        return (md5gen<item_t > (item) == hash);
    }

    extern std::string md5gen(const std::string & data);
    extern std::string md5gen(const fs::path & file);
};



#endif	/* _SERIALHASH_H */


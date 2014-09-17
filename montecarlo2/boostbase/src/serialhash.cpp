#include "serialhash.h"

namespace boostbase
{
    std::string md5gen(const std::string & data){
        hashwrapper * md5 = new md5wrapper();
        std::string hash = md5->getHashFromString(data);
        delete md5;
        return hash;
    }
    std::string md5gen(const fs::path & file){
        hashwrapper * md5 = new md5wrapper();
        std::string hash = md5->getHashFromFile(file.string());
        delete md5;
        return hash;
    }
}

/* 
 * File:   saveserializable.h
 * Author: karol
 *
 * Created on 3 listopad 2009, 19:49
 */

#ifndef _SERIALSAVELOAD_H
#define	_SERIALSAVELOAD_H

#include "serializer.h"
#include "serialhash.h"
#include "boost.h"
#include "exceptions.h"
#include "stringsaveload.h"

namespace boostbase {

    /*
     * zapisywanie obiektu do pliku i czytanie
     */
    template<class item_t>
    void serialsave(const item_t & item, const fs::path & file) {
        if (!file.parent_path().empty() && !fs::exists(file.parent_path()))
            fs::create_directories(file.parent_path());
        if (fs::exists(file.parent_path()) && !fs::is_directory(file.parent_path()))
            throw exception::file_is_directory();
        std::ofstream o_file;
        o_file.exceptions(std::ofstream::badbit | std::ofstream::failbit);
        try {
            o_file.open(file.string().c_str(), std::ios_base::out | std::ios_base::binary);
        } catch (std::ofstream::failure e) {
            if (o_file.fail())
                throw exception::failed_opening_file();
            if (o_file.bad())
                throw exception::failed_opening_file();
        }

        std::string out = boostbase::serialsave<item_t>(item);
        o_file << out;
        o_file.close();
    }

    template<class item_t>
    item_t serialload(const fs::path & file,const std::string & filehash="") {
        if (!fs::exists(file))
            throw exception::file_not_found();
        std::ifstream i_file;
        
        i_file.exceptions(std::ifstream::badbit);
        try {
            i_file.open(file.string().c_str(), std::ios_base::in | std::ios_base::binary);
        } catch (std::ifstream::failure e) {
            if (i_file.eof())
                throw exception::unexpected_eof();
            if (i_file.bad() || i_file.fail())
                throw exception::failed_opening_file();
        }
        
        //i_file.open(file.string().c_str(), std::ios_base::in | std::ios_base::binary);
        if(filehash!=""){
            if(md5gen(file)!=filehash){
                std::cout << md5gen(file) << " != " << filehash << std::endl;
                throw exception::wrong_md5();
            }
        }

        std::string data;
        try {
            io::copy(i_file,io::back_inserter(data));
            //i_file >> data;
        } catch (std::ofstream::failure e) {
            throw exception::failed_reading_file();
        }
        item_t item = boostbase::serialload<item_t>(data);
        try {
            i_file.close();
        } catch (std::ofstream::failure e) {
            throw exception::failed_reading_file();
        }

        return item;


    }

}


#endif	/* _SERIALSAVELOAD_H */


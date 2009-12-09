/* 
 * File:   stringsaveload.h
 * Author: karol
 *
 * Created on 8 listopad 2009, 14:26
 */

#ifndef _STRINGSAVELOAD_H
#define	_STRINGSAVELOAD_H

#include "serializer.h"

namespace boostbase {

    /**
     * zamiana obiektu na string i odwrotnie
     * 
     * jakby były wątpliwości, to jest tutaj kompresja i dekompresja przez zlib
     */
    template<class item_t>
    std::string serialsave(item_t & item) {
        std::string output;
        io::filtering_ostream oo;
        oo.push(io::zlib_compressor(io::zlib_params(9,io::zlib::deflated,15,9,io::zlib::default_strategy,false)));
        oo.push(io::back_inserter(output));
        //oo.push(output);
        ofilterserializer of(oo);
        of | item;
        return output;
    }

    template<class item_t>
    item_t serialload(const std::string & s) {
        io::filtering_istream ii;
        ii.push(io::zlib_decompressor());
        ii.push(boost::make_iterator_range(s));
        //ii.push(s);
        ifilterserializer ifser(ii);
        item_t item;
        ifser | item;
        return item;
    }

};


#endif	/* _STRINGSAVELOAD_H */


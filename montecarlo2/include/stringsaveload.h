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

    /// serialization
    template<class item_t>
    int size(const item_t & item) {
        io::filtering_ostream oo;
        oo.push(io::counter());
        oo.push(io::null_sink());
        ofilterserializer of(oo);
        of | const_cast< typename remove_const<item_t>::type &>(item);
        return oo.component<0, io::counter>()->characters();
    }

    template<class item_t>
    std::string serialsave(const item_t & item) {
        int i_size = size<item_t>(item);

        std::string output;
        io::filtering_ostream oo;
        oo.push(io::zlib_compressor(io::zlib_params(9,io::zlib::deflated,15,9,io::zlib::default_strategy,false),2*i_size),3*i_size);
        oo.push(io::back_inserter(output));
        //oo.push(output);
        ofilterserializer of(oo);
        of | const_cast< typename remove_const<item_t>::type &>(item);
        return output;
    }

    template<class item_t>
    item_t serialload(const std::string & s) {
        int b_size = s.size();

        io::filtering_istream ii;
        ii.push(io::zlib_decompressor(io::zlib_params(9,io::zlib::deflated,15,9,io::zlib::default_strategy,false),2*b_size),3*b_size);
        ii.push(boost::make_iterator_range(s));
        //ii.push(s);
        ifilterserializer ifser(ii);
        item_t item;
        ifser | item;
        return item;
    }

}


#endif	/* _STRINGSAVELOAD_H */


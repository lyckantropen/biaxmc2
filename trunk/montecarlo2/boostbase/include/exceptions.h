/* 
 * File:   bb_exceptions.h
 * Author: karol
 *
 * Created on 3 listopad 2009, 21:11
 */

#ifndef _BB_EXCEPTIONS_H
#define	_BB_EXCEPTIONS_H

#include "std.h"

namespace boostbase{ namespace exception {

    class file_not_found:public std::exception {};
    class file_is_directory:public std::exception {};
    class failed_reading_file:public std::exception {};
    class unexpected_eof:public std::exception {};
    class failed_opening_file:public std::exception {};
    class zlib_stream_read_error:public std::exception {};
    class zlib_courrupted_input_stream:public std::exception {};
    class zlib_buffer_read_error:public std::exception {};
    class zlib_stream_write_error:public std::exception {};
    class zlib_courrupted_output_stream:public std::exception {};
    class zlib_buffer_write_error:public std::exception {};
    class sqlite_failed_opening_db:public std::exception {};
    class sqlite_failed_creating_db:public std::exception {};
    class sqlite_command_error:public std::exception {};
    class wrong_md5:public std::exception {};

};};


#endif	/* _BB_EXCEPTIONS_H */


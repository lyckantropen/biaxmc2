/* 
 * File:   boost_headers.h
 * Author: karol
 *
 * Created on 29 październik 2009, 16:30
 */

#ifndef _BOOST_HEADERS_H
#define	_BOOST_HEADERS_H

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filter/counter.hpp>
#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>

#define foreach BOOST_FOREACH

namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace pt = boost::posix_time;
namespace gr = boost::gregorian;
namespace po = boost::program_options;
namespace assign = boost::assign;
using boost::shared_ptr;
using boost::make_shared;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;


#endif	/* _BOOST_HEADERS_H */


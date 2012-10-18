/* 
 * File:   serializer.h
 * Author: karol
 *
 * Created on 3 listopad 2009, 16:53
 */

#ifndef _I_SERIALIZER_H
#define	_I_SERIALIZER_H

#include "std.h"
#include "boost.h"
//#include "zfile.h"

namespace boostbase {

    class serializer {
    public:
        virtual void operator|(char &) = 0;
        virtual void operator|(double &) = 0;
        virtual void operator|(float &) = 0;
        virtual void operator|(int &) = 0;
        virtual void operator|(unsigned int &) = 0;
        virtual void operator|(long int &) = 0;
        virtual void operator|(short int &) = 0;
        virtual void operator|(bool &) = 0;
        virtual void operator|(std::string &) = 0;
        virtual void operator|(std::list<float> &) = 0;
        virtual void operator|(std::list<double> &) = 0;
        virtual void operator|(std::list<std::string> &) = 0;
        virtual void operator|(std::vector<float> &) = 0;
        virtual void operator|(std::vector<double> &) = 0;
        virtual void operator|(std::vector<int> &) = 0;
        virtual void operator|(std::vector<std::string> &) = 0;
        virtual void operator|(std::valarray<double> &) = 0;
        virtual void flush() = 0;
    };

    template<class stream_t>
    class outserializer : public serializer {
        stream_t & s;
    public:

        outserializer(stream_t & d) : s(d) {
        }

        outserializer(const outserializer & c) : s(c.s) {
        }

        virtual void operator|(char & c) {
            s.write((char*) & c, sizeof (char));
        }

        virtual void operator|(float & c) {
            s.write((char*) & c, sizeof (float));
        }

        virtual void operator|(double & c) {
            s.write((char*) & c, sizeof (double));
        }

        virtual void operator|(int & c) {
            s.write((char*) & c, sizeof (int));
        }

        virtual void operator|(unsigned int & c) {
            s.write((char*) & c, sizeof (unsigned int));
        }

        virtual void operator|(short int & c) {
            s.write((char*) & c, sizeof (short int));
        }
        
        virtual void operator|(bool & c) {
            s.write((char*) & c, sizeof (bool));
        }


        virtual void operator|(long int & c) {
            s.write((char*) & c, sizeof (long int));
        }

        virtual void operator|(std::string & c) {
            unsigned int len = c.length();
            operator|(len);
            for (int i = 0; i < len; i++) {
                operator|(c[i]);
            }
        }

        template<class t>
        void operator|(std::vector<t> & v) {
            unsigned int size = v.size();
            operator|(size);

            foreach(t & item, v) {

                operator|(item);
            }
        }
        template<class t>
        void operator|(std::valarray<t> & v) {
            unsigned int size = v.size();
            operator|(size);
            //std::cout << "out: " << size << std::endl;
            for(int i=0;i<size;i++) {
                operator|(v[i]);
            }
        }

        template<class t>
        void operator|(std::list<t> & l) {
            unsigned int size = l.size();
            operator|(size);

            foreach(t & item, l) {

                operator|(item);
            }
        }

        virtual void operator|(std::list<float> & c) {
            operator|<float>(c);
        }

        virtual void operator|(std::list<double> & c) {
            operator|<double>(c);
        }

        virtual void operator|(std::list<std::string> & c) {
            operator|<std::string > (c);
        }

        virtual void operator|(std::vector<float> & c) {
            operator|<float>(c);
        }

        virtual void operator|(std::vector<double> & c) {
            operator|<double>(c);
        }
        virtual void operator|(std::vector<int> & c) {
            operator|<int>(c);
        }

        virtual void operator|(std::vector<std::string> & c) {
            operator|<std::string > (c);
        }
        virtual void operator|(std::valarray<double> & c) {
            operator|<double> (c);
        }

        virtual void flush() {
            s.flush();
        }
    };

    template<class stream_t>
    class inserializer : public serializer {
        stream_t & s;
    public:

        inserializer(stream_t & d) : s(d) {
        }

        inserializer(const inserializer & c) : s(c.s) {
        }

        virtual void operator|(char & c) {
            s.read((char*) & c, sizeof (char));
        }

        virtual void operator|(float & c) {
            s.read((char*) & c, sizeof (float));
        }

        virtual void operator|(double & c) {
            s.read((char*) & c, sizeof (double));
        }

        virtual void operator|(int & c) {
            s.read((char*) & c, sizeof (int));
        }

        virtual void operator|(unsigned int & c) {
            s.read((char*) & c, sizeof (unsigned int));
        }

        virtual void operator|(short int & c) {
            s.read((char*) & c, sizeof (short int));
        }
        
        virtual void operator|(bool & c) {
            s.read((char*) & c, sizeof (bool));
        }

        virtual void operator|(long int & c) {
            s.read((char*) & c, sizeof (long int));
        }

        virtual void operator|(std::string & c) {
            unsigned int length;
            operator|(length);
            c.clear();
            for (int i = 0; i < length; i++) {
                char C;
                s.read(&C, sizeof (char));
                c.push_back(C);
            }
            c.push_back('\0');
        }

        template<class t>
        void operator|(std::vector<t> & v) {
            unsigned int size;
            operator|(size);
            v.reserve(size);
            for (int i = 0; i < size; i++) {
                t item;
                operator|(item);
                v.push_back(item);
            }
        }
        template<class t>
        void operator|(std::valarray<t> & v) {
            unsigned int size;
            operator|(size);
            v.resize(size);
            //std::cout << "in: "<< size << std::endl;
            for (int i = 0; i < size; i++) {
                t item;
                operator|(item);
                v[i]=item;
            }
        }
        template<class t>
        void operator|(std::list<t> & l) {
            unsigned int size;
            operator|(size);
            for (int i = 0; i < size; i++) {
                t item;
                operator|(item);
                l.push_back(item);
            }
        }

        virtual void operator|(std::list<float> & c) {
            operator|<float>(c);
        }

        virtual void operator|(std::list<double> & c) {
            operator|<double>(c);
        }

        virtual void operator|(std::list<std::string> & c) {
            operator|<std::string > (c);
        }

        virtual void operator|(std::vector<float> & c) {
            operator|<float>(c);
        }

        virtual void operator|(std::vector<double> & c) {
            operator|<double>(c);
        }
        virtual void operator|(std::vector<int> & c) {
            operator|<int>(c);
        }
        virtual void operator|(std::valarray<double> & c) {
            operator|<double>(c);
        }

        virtual void operator|(std::vector<std::string> & c) {
            operator|<std::string > (c);
        }

        virtual void flush() {
        }
    };
    
    typedef outserializer<std::stringstream> osstreamserializer;
    typedef inserializer<std::stringstream> isstreamserializer;
    typedef outserializer<std::ofstream> ofileserializer;
    typedef inserializer<std::ifstream> ifileserializer;
    typedef outserializer<io::filtering_ostream> ofilterserializer;
    typedef inserializer<io::filtering_istream> ifilterserializer;

};

#endif	/* _I_SERIALIZER_H */


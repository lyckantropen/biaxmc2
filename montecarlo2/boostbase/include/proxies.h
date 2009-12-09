/* 
 * File:   proxies.h
 * Author: karol
 *
 * Created on 8 listopad 2009, 19:28
 */

#ifndef _PROXIES_H
#define	_PROXIES_H
namespace boostbase {

    typedef std::pair<std::string, std::string> pair_t;

    ///predykat do sortowania wektora pair_t
    bool pair_t_pred(const pair_t & a, const pair_t & b) {
        return a.first <= b.first;
    }
    
    /**
     * obiekt proxy pozwalający na generowanie listy warunków za pomocą operatora ()
     */
    class pair_t_proxy {
        std::vector<pair_t> values;
    public:
        pair_t_proxy() {}
        pair_t_proxy(const std::vector<pair_t> & val):values(val) {}
        template<class value_type>
        pair_t_proxy operator()(const std::string & field, const value_type & val){
            std::stringstream s;
            s << typeid(val).name() << ":" << val;
            values.push_back(pair_t(field,s.str()));
            return pair_t_proxy(values);
        }
        operator const std::vector<pair_t> & () const {
            return values;
        }
    };

    ///tworzenie listy warunków: where("wielkość1","wartość1")("wielkość2","wartość2")("wielkość3","wartość3")..
    template<class value_type>
    pair_t_proxy where(const std::string & field, const value_type & val) {
        std::vector<pair_t> values;
        std::stringstream s;
        s << typeid(val).name() << ":" << val;
        values.push_back(pair_t(field,s.str()));
        return pair_t_proxy(values);
    }
    template pair_t_proxy where<double>(const std::string &, const double &);
    template pair_t_proxy where<float>(const std::string &, const float &);
    template pair_t_proxy where<int>(const std::string &, const int &);
    template pair_t_proxy where<std::string>(const std::string &, const std::string &);
    template pair_t_proxy pair_t_proxy::operator()<double>(const std::string &,const double &);
    template pair_t_proxy pair_t_proxy::operator()<float>(const std::string &, const float &);
    template pair_t_proxy pair_t_proxy::operator()<int>(const std::string &, const int &);
    template pair_t_proxy pair_t_proxy::operator()<std::string>(const std::string &, const std::string &);

    struct tween_t {
        std::string field,begin,end;
        tween_t(const std::string & a, const std::string & b, const std::string & c):
        field(a),begin(b),end(c)
        {}
        bool operator!=(const tween_t & c) const {
            return (field!=c.field && begin!=c.begin && end!=c.end);
        }
    };

    class tween_t_proxy {
        std::vector<tween_t> values;
    public:
        tween_t_proxy() {}
        tween_t_proxy(const std::vector<tween_t> & val):values(val) {}
        template<class value_type>
        tween_t_proxy operator()(const std::string & field, const value_type & begin,const value_type & end ){
            std::stringstream s1,s2;
            s1 << begin;
            s2 << end;
            values.push_back(tween_t(field,s1.str(),s2.str()));
            return tween_t_proxy(values);
        }
        operator const std::vector<tween_t> & () const {
            return values;
        }
    };

    template<class value_type>
    tween_t_proxy between(const std::string & field, const value_type & begin, const value_type & end){
        std::vector<tween_t> values;
        std::stringstream s1,s2;
        s1 << begin;
        s2 << end;
        values.push_back(tween_t(field,s1.str(),s2.str()));
        return tween_t_proxy(values);
    }
    template tween_t_proxy between<double>(const std::string &, const double &, const double &);
    template tween_t_proxy between<float>(const std::string &, const float &, const float &);
    template tween_t_proxy between<int>(const std::string &, const int &, const int &);
    template tween_t_proxy between<std::string>(const std::string &, const std::string &, const std::string &);
    template tween_t_proxy tween_t_proxy::operator()<double>(const std::string &, const double &, const double &);
    template tween_t_proxy tween_t_proxy::operator()<float>(const std::string &, const float &, const float &);
    template tween_t_proxy tween_t_proxy::operator()<int>(const std::string &, const int &, const int &);
    template tween_t_proxy tween_t_proxy::operator()<std::string>(const std::string &, const std::string &, const std::string &);

};


#endif	/* _PROXIES_H */


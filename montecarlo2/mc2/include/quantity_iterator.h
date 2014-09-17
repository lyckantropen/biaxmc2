/*
 * File:   valuearray.h
 * Author: karol
 *
 * Created on 9 listopad 2009, 18:12
 */

#ifndef _VALUEARRAY_H
#define _VALUEARRAY_H

#include "std.h"
#include "boost.h"


/*
 * klasa quantity_iterator zastępuje klasę Krata
 * niestety nie jest ona thread-safe
 */

template<class value_type>
class ranged_quantity
{
    value_type & value;
    value_type low;
    value_type high;
    value_type interval;
public:
    typedef enum { up, down } dir_t;
private:
    dir_t direction;
public:
    ranged_quantity(value_type & val, const value_type & l, const value_type & h, const value_type & i):
        value(val), low(l), high(h), interval(i)
    {
        if(low > high)
            direction = down;
        if(high > low)
            direction = up;
        value = low;
    }
    ranged_quantity(const ranged_quantity & source): value(source.value)
    {
        low = source.low;
        high = source.high;
        interval = source.interval;
        if(low > high)
            direction = down;
        if(high > low)
            direction = up;
    }
    const ranged_quantity & operator=(const ranged_quantity & source)
    {
        value = source.value;
        low = source.low;
        high = source.high;
        interval = source.interval;
        if(low > high)
            direction = down;
        if(high > low)
            direction = up;
        return *this;
    }

    operator value_type()
    {
        return value;
    }
    bool operator++(int)
    {
        if(finished())
            return false;
        if(direction == up)
            value += interval;
        if(direction == down)
            value -= interval;
        return true;
    }
    bool finished()
    {
        if(direction == up)
            return (value + interval) > end();
        else return (value - interval) < end();
    }
    void reset()
    {
        value = low;
    }
private:
    value_type begin() const
    {
        if(direction == up)
            return low;
        else return high;
    }
    value_type end() const
    {
        if(direction == up)
            return high - interval;
        else return high + interval;
    }
};

template<class quantity_t, class value_type>
class quantity_iterator
{
    template<class v> friend class value_array;
    typedef typename std::list<quantity_t>::iterator     iterator;
    std::list<quantity_t>                       quantities;
    iterator                                current;
public:
    quantity_iterator(const std::list<quantity_t> & q):
        quantities(q)
    {
        // z niewiadomego powodu to tutaj musi być
        current = quantities.end();
    }
    void add(const quantity_t & q)
    {
        quantities.push_back(q);
    }
    bool operator++(int)
    {
        iterator next = quantities.end();
        next = find_last_unfinished(next);
        if(next == quantities.begin() && next->finished())
            return false;
        (*next)++;
        current = next;
        next++;
        for(iterator i = next; i != quantities.end(); i++)
            i->reset();
        if(next == quantities.end())
            next--;
        return true;
    }
    operator value_type() const
    {
        return *current;
    }
private:
    iterator find_last_unfinished(iterator & it)
    {
        if(!it->finished() || it == quantities.begin())
            return it;
        else
        {
            it--;
            return find_last_unfinished(it);
        }
    }
};

typedef ranged_quantity<double> ranged_double;
typedef ranged_quantity<float> ranged_float;
typedef ranged_quantity<int> ranged_int;

typedef quantity_iterator<ranged_quantity<double>, double>   double_iterator;
typedef quantity_iterator<ranged_quantity<float>, float>    float_iterator;
typedef quantity_iterator<ranged_quantity<int>, int>      int_iterator;

template<class value_type>
class ranged_quantity_proxy
{
    std::list<ranged_quantity<value_type> > values;
public:
    ranged_quantity_proxy(const std::list<ranged_quantity<value_type> > & v): values(v) {}
    ranged_quantity_proxy operator()(value_type & val, const value_type & a, const value_type & b, const value_type & c)
    {
        values.push_back(ranged_quantity<value_type>(val, a, b, c));
        return ranged_quantity_proxy(values);
    }
    operator std::list<ranged_quantity<value_type> > () const
    {
        return values;
    }
};

template<class value_type>
ranged_quantity_proxy<value_type> with(value_type & val, const value_type & a, const value_type & b, const value_type & c)
{
    std::list<ranged_quantity<value_type> > values;
    values.push_back(ranged_quantity<value_type>(val, a, b, c));
    return ranged_quantity_proxy<value_type>(values);
};
template ranged_quantity_proxy<double> with<double>(double &, const double &, const double &, const double &);



#endif  /* _VALUEARRAY_H */


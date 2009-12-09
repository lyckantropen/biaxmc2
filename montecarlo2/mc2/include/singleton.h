/* 
 * File:   singleton.h
 * Author: karol
 *
 * Created on 17 listopad 2009, 12:58
 */

#ifndef _SINGLETON_H
#define	_SINGLETON_H

template<class t>
class Singleton {
    static  t   *   instance;
protected:
    Singleton(){}
public:
    static t * Instance() {
        if(instance==NULL){
            instance=new t;
        }
        return instance;
    }
};

#endif	/* _SINGLETON_H */


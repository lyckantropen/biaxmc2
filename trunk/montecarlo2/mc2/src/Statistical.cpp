#include "Statistical.h"

Value CalculateFluctuation(const vect & variable,const int & acc_idx){
    int size=acc_idx+1;
    if(acc_idx==0)
        size=variable.size();

    vect fluct3(size);
    for(int i=0;i<(size);i++){
        vect e(size);
        vect e2(size);
        for(int j=0;j<(size);j++){
            int t = (size)*random01();
            e[j]=variable[t];
            e2[j]=variable[t]*variable[t];
        }
        fluct3[i]=(double(Mean(e2))-double(Mean(e))*double(Mean(e)));
    }
    return Mean(fluct3);
}

///
/// Value implementation
///

Value::Value():val(0),err(0) {}

Value::Value(const double &v, const double &e){
    val=v;
    err=e;
}

Value::Value(const double &v){
    val=v;
    err=0.0;
}

const double &Value::Error() const {
    return err;
}

const Value Value::operator=(const double &v){
    val=v;
    return *this;
}

const Value Value::operator+(const Value &v) const {
    return Value(val+v.val,err+v.err);
}

const Value Value::operator-(const Value &v) const {
    return Value(val-v.val,err+v.err);
}

const Value Value::operator*(const Value &v) const {
    return Value(val*v.val,std::abs(err*v.val)+std::abs(v.err*val));
}

const Value Value::operator/(const Value &v) const {
    return Value(val/v.val,std::abs(err/v.val)+std::abs(v.err/val/val));
}

bool Value::operator==(const Value &v){
    return val==v.val;
}

const std::string Value::Print() const {
    std::stringstream s;
    s << std::setprecision(2);
    s << val << "("<<err<<","<< std::abs(err/val*100) << "%)";
    return s.str();
}

const std::string Value::TableForm() const {
    std::stringstream s;
    s << val << "\t" << err ;
    return s.str();
}

const std::string Value::MathematicaForm() const {
    std::stringstream s;
    s << std::setprecision(15) << std::fixed;
    s << "{" << val << ","<<err<<"}";
    return s.str();
}

bool Value::operator==(const double &v){
    return val==v;
}

Value::operator const double &() const {
    return val;
}

Value::operator double &(){
    return val;
}

std::ostream & operator<<(std::ostream & s, Value & v){
    s << double(v) << " " << v.Error() ;
    return s;
}

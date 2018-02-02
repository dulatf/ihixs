#ifndef RESULT_PAIR_H
#define RESULT_PAIR_H

#include "math.h"
using namespace std;
class ResultPair{
public:
    /// Empty Constructor
    ResultPair():_val(0.0),_err(0.0){};
    /// Constructor
    ResultPair(const double& val, const double& err):_val(val),_err(err){};
    
    /// Copy constructor
    ResultPair(const ResultPair& that) :
    _val(that.val()),_err(that.err())
    {}
    
    /// Move constructor
    ResultPair(ResultPair&& that) :
    _val(that.val()),_err(that.err())
    {}
    
    /// Destructor
    ~ResultPair()
    {}
    
    /// @}
    
    /// \name Input/output functions
    /// @{
    
    /// Assignment operator
    ResultPair& operator=(const ResultPair& that)
    {
        _val = that.val();
        _err = that.err();
        return (*this);
    }
    
    /// Move assignment operator
    ResultPair& operator=(ResultPair&& that)
    {
        _val = that.val();
        _err = that.err();
        return (*this);
    }
    
    /// Return the value
    double val() const {return _val;}
    
    /// Return the error
    double err() const {return _err;}
    
    
    /// @}
    
    /// \name Operations
    /// @{
    
    /// Multiply this ResultPair with a number
    ResultPair operator*(const double& that) const
    {
        return ResultPair(that*_val,that*_err);
    }
    
    /// Multiply the ResultPair v with a number k
    friend ResultPair operator*(const double k, const ResultPair& v)
    {
        return ResultPair(k*v.val(),k*v.err());
    }
    
    /// Multiply this ResultPair with another Resultpair
    ResultPair operator*(const ResultPair& that) const
    {
        return ResultPair(that.val()*_val,
                          sqrt(pow(that.val()*_err,2.)+pow(_val*that.err(),2.))
                          );
    }
    
    /// Divide this ResultPair by another Resultpair
    ResultPair operator/(const ResultPair& that) const
    {
        if ( _val == 0.0 ) return ResultPair(0.0,0.0);
        else
            return ResultPair(_val / that.val(),
                          fabs(_val / that.val()) *
                          sqrt(pow(that.err()/that.val(),2.)+pow(_err/_val,2.))
                          );
    }
    
    /// Add this ResultPair to another one
    ResultPair operator+(const ResultPair& that) const
    {
        return ResultPair(_val+that.val(),sqrt(pow(_err,2.)+pow(that.err(),2.)));
    }
    
    /// Subtract from this ResultPair  another one
    ResultPair operator-(const ResultPair& that) const
    {
        return ResultPair(_val-that.val(),sqrt(pow(_err,2.)+pow(that.err(),2.)));
    }
    
    void TakeAbsoluteValue() {if (_val < 0) _val = -_val;}
    
    friend ostream& operator<<(ostream& stream, const ResultPair& );
    
private:
    double _val;
    double _err;
};

#endif

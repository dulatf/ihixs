

#ifndef AS_SERIES_H
#define AS_SERIES_H

#include "string"
#include "vector"
#include "iostream"
#include "iomanip"
using namespace std;

#include "result_pair.h"

/**
 *
 * \class AsSeries
 * \brief A representation of an expansion in a_s/pi
 * \brief a_s^k * c[k] + a_s^(k+1) * c[k+1] + ... + a_s^[k+m] c[k+m]
 * \brief Each c[k] is a ResultPair, i.e. a double value with an error
 * \brief multiplication, addition and truncation are defined
 *
 */

class AsSeries{
public:
    /// Constructor of empty series
    AsSeries() :
    _starting_exponent(0)
    {_terms.push_back(ResultPair(0.0,0.0));}
    /// Constructor with one term
    AsSeries(int starting_exponent,const double& val1);
    /// Constructor with one term
    AsSeries(int starting_exponent,const ResultPair& rp1);
    /// Constructor with two terms
    AsSeries(int starting_exponent  ,const double& val1
                                    ,const double& val2);
    /// Constructor with three terms
    AsSeries(int starting_exponent,const double& val1
                                    ,const double& val2
                                    ,const double& val3);
    /// Constructor with four terms
    AsSeries(int starting_exponent,const double& val1
             ,const double& val2
             ,const double& val3
             ,const double& val4);
    
    /// Constructor from vector<ResultPairs>
    AsSeries(int start_exp,const vector<ResultPair>& terms) :
    _terms(terms),_starting_exponent(start_exp)
    {}
    
    /// Copy constructor
    AsSeries(const AsSeries& that) :
    _terms(that.terms()), _starting_exponent(that.starting_exponent())
    {
    }
    
    /// Move constructor
    AsSeries(AsSeries&& that) :
    _terms(that.terms()), _starting_exponent(that.starting_exponent())
    {
    }
    
    /// Destructor
    ~AsSeries()
    {}
    
    /// @}
    
    /// \name Assignment operators
    /// @{
    /// Assignment operator
    AsSeries& operator=(const AsSeries& that)
    {
        _terms = that.terms();
        _starting_exponent=that.starting_exponent();
        return (*this);
    }
    
    /// Move assignment operator
    AsSeries& operator=(AsSeries&& that)
    {
        _terms = that.terms();
        _starting_exponent=that.starting_exponent();
        return (*this);
    }
    
    /// @}
    
    /// \name Input/output functions
    /// @{
    
    /// Return starting power of a_s
    int starting_exponent()const
    {return _starting_exponent;}
    /// Return ending power of a_s
    int ending_exponent() const
    {return int(_terms.size())+_starting_exponent-1;}
    
    /// Return the term of order n
    ResultPair term_of_order(int n) const
    {
    if (n>_starting_exponent-1 and n<ending_exponent()+1)
        return _terms[n-_starting_exponent];
    else
        return ResultPair(0.0,0.0);
    }
    
    /// Return the terms as a vector<ResultPairs>
    vector<ResultPair> terms() const {return _terms;}
    
    ///Return the sum of terms
    ResultPair AddUp(){ResultPair res(0.0,0.0);for (int i=0;i<_terms.size();i++) res = res + _terms[i]; return res;}
    
    /// @}
    
    /// \name Operations
    /// @{
    
    /// Multiply this series with a number
    AsSeries operator*(const double& that) const;

    /// Multiply this series with a ResultPair
    AsSeries operator*(const ResultPair& that) const;
    
    /// Multiply this series with a ResultPair
    friend AsSeries operator*(const ResultPair& that, const AsSeries& theother);
    
    /// Multiply each term with a_s/pi to the proper power
    void MultiplyAs(const double& as_pi);
    
    
    /// Shift the starting exponent by n (equivalent to multiplying by a_s to n)
    void ShiftStartingExponentBy(int n){_starting_exponent = _starting_exponent+n;}
    
    /// Multiply each term by a_s^n where n is the exponent of the term
    /// and a_s is itself a series (effectively, set a_s to be a series in a_s)
    /// used to evolve as(muf)-> as(mur)
    AsSeries MultiplyBySeries(const AsSeries& as_mur);

    /// Return this AsSeries to the power k
    AsSeries ToThe(int k) const;
    
    /// Multiply the AsSeries v with a number k
    friend AsSeries operator*(const double k, const AsSeries& v)
    {
        return v.operator*(k);
    }
    
    /// Add this AsSeries to another one
    AsSeries operator+(const AsSeries& that) const;
    
    /// Multiply two AsSeries
    AsSeries operator*(const AsSeries& that) const;
    
    /// Truncate the series to order a_s^n (terms with k>n are discarded)
    void Truncate(int n);
    
    /// Is the series trivially zero (i.e. is it empty)?
    bool IsZero(){return _terms.empty();}
    ///@}
    /// @{
    /// \name Printing
    friend ostream& operator<<(ostream& stream, const AsSeries& st);

    string PrintCumulative();
    
    //@}
private:
    vector<ResultPair> _terms;
    int _starting_exponent;
};

#endif

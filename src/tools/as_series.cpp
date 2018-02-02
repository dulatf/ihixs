#include "as_series.h"
#include "iomanip"
#include<sstream>
using namespace std;



AsSeries::AsSeries(int starting_exponent,const double& val1)
{
    _starting_exponent = starting_exponent;
    _terms.push_back(ResultPair(val1,0.0));
}

AsSeries::AsSeries(int starting_exponent,const ResultPair& val1)
{
    _starting_exponent = starting_exponent;
    _terms.push_back(val1);
}

AsSeries::AsSeries(int starting_exponent,const double& val1, const double&val2)
{
    _starting_exponent = starting_exponent;
    _terms.push_back(ResultPair(val1,0.0));
    _terms.push_back(ResultPair(val2,0.0));
}

AsSeries::AsSeries(int starting_exponent,const double& val1, const double&val2, const double& val3)
{
    _starting_exponent = starting_exponent;
    _terms.push_back(ResultPair(val1,0.0));
    _terms.push_back(ResultPair(val2,0.0));
    _terms.push_back(ResultPair(val3,0.0));
}

AsSeries::AsSeries(int starting_exponent,const double& val1, const double&val2,const double& val3, const double& val4)
{
    _starting_exponent = starting_exponent;
    _terms.push_back(ResultPair(val1,0.0));
    _terms.push_back(ResultPair(val2,0.0));
    _terms.push_back(ResultPair(val3,0.0));
    _terms.push_back(ResultPair(val4,0.0));
}

AsSeries AsSeries::operator*(const double& that) const
{
    vector<ResultPair> terms = _terms;
    for (int i=0;i<terms.size();i++) terms[i] = that*terms[i];
    return AsSeries(_starting_exponent, terms);
}

AsSeries AsSeries::operator*(const ResultPair& that) const
{
    vector<ResultPair> terms = _terms;
    for (int i=0;i<terms.size();i++)
    {
        terms[i] = that*terms[i];
    }
    return AsSeries(_starting_exponent, terms);
}

void AsSeries::MultiplyAs(const double& as_pi)
{
    //cout<<__func__<<endl;
    for (int k=_starting_exponent;k<ending_exponent()+1;k++)
    {
        _terms[k-_starting_exponent] = _terms[k-_starting_exponent]*pow(as_pi,double(k));
    }
}

AsSeries AsSeries::MultiplyBySeries(const AsSeries& as_mur)
{
    AsSeries  res;
    for (int k=_starting_exponent;k<ending_exponent()+1;k++)
    {
        res = res + _terms[k-_starting_exponent] * as_mur.ToThe(k);
    }
    return res;
}

AsSeries operator*(const ResultPair& that, const AsSeries& theother) 
{
    return theother*that;
}



AsSeries AsSeries::ToThe(int k) const
{
    AsSeries res =  *this;
    for (int i=0;i<k-1;i++)
    {
        res = res * (*this);
    }
    return res;
}


AsSeries AsSeries::operator+(const AsSeries& that) const
{
    int k_initial = min(_starting_exponent,that.starting_exponent());
    int k_final = max(ending_exponent(),that.ending_exponent());
    vector<ResultPair> terms;
    
    for (int k=k_initial;k<k_final+1;k++)
    {
        terms.push_back(term_of_order(k)+that.term_of_order(k));
    }
    return AsSeries(k_initial, terms);
}

AsSeries AsSeries::operator*(const AsSeries& that) const
{
    int k_initial = _starting_exponent+that.starting_exponent();
    int k_final = ending_exponent()+that.ending_exponent();
    vector<ResultPair> terms;
    
    for (int k=k_initial;k<k_final+1;k++)
    {
        ResultPair newterm(0.0,0.0);
        for (int m=_starting_exponent;m<ending_exponent()+1;m++)
        {
            newterm = newterm + term_of_order(m) * that.term_of_order(k-m);
        }
        
        terms.push_back(newterm);
    }
    return AsSeries(k_initial, terms);
}

void AsSeries::Truncate(int n)
{
    vector<ResultPair> newterms;
    for (int i=0;i<_terms.size();i++)
    {
        if (i+_starting_exponent<n+1) newterms.push_back(_terms[i]);
    }
    _terms = newterms;
}

ostream& operator<<(ostream& stream, const AsSeries& st)
{
    for (int i=st.starting_exponent();i<st.ending_exponent()+1;i++)
    {
        if (fabs(st.term_of_order(i).val())>1e-15)
            stream<<"a^"<<i<<":"<<st.term_of_order(i);
        
    }
    return stream;
}


string AsSeries::PrintCumulative()
{


    stringstream s(stringstream::out);
    ResultPair res;
    for (int i=starting_exponent();i<ending_exponent()+1;i++)
    {
        
        if (fabs(term_of_order(i).val())>1e-15){
            res = res+term_of_order(i);
            s<<" "<<res;
        }
    }
    return s.str();
}
















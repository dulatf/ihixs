#ifndef CPP_CHAPLIN_H
#define CPP_CHAPLIN_H

#include <complex>
#include <vector>
using namespace std;


namespace chaplin{
    complex<double> HPL(int n1, const complex<double>& x);
    complex<double> HPL(int n1, int n2, const complex<double>& x);
    complex<double> HPL(int n1, int n2, int n3, const complex<double>& x);
    complex<double> HPL(int n1, int n2, int n3, int n4, const complex<double>& x);
    
    complex<double> basis2x1(const complex<double>& x);
    complex<double> basis2x2(const complex<double>& x);
    complex<double> basis2x3(const complex<double>& x);

    complex<double> basis3x1(const complex<double>& x);
    complex<double> basis3x2(const complex<double>& x);
    complex<double> basis3x3(const complex<double>& x);
    complex<double> basis3x4(const complex<double>& x);
    complex<double> basis3x5(const complex<double>& x);
    complex<double> basis3x6(const complex<double>& x);
    complex<double> basis3x7(const complex<double>& x);
    complex<double> basis3x8(const complex<double>& x);
    
    complex<double> basis1(const complex<double>& x);
    complex<double> basis2(const complex<double>& x);
    complex<double> basis3(const complex<double>& x);
    complex<double> basis4(const complex<double>& x);
    complex<double> basis5(const complex<double>& x);
    complex<double> basis6(const complex<double>& x);
    complex<double> basis7(const complex<double>& x);
    complex<double> basis8(const complex<double>& x);
    complex<double> basis9(const complex<double>& x);
    complex<double> basis10(const complex<double>& x);
    complex<double> basis11(const complex<double>& x);
    complex<double> basis12(const complex<double>& x);
    complex<double> basis13(const complex<double>& x);
    complex<double> basis14(const complex<double>& x);
    complex<double> basis15(const complex<double>& x);
    complex<double> basis16(const complex<double>& x);
    complex<double> basis17(const complex<double>& x);
    complex<double> basis18(const complex<double>& x);

    complex<double> bsh21m1_outside_1(const complex<double>& x);
    complex<double> bsh21m1_outside_2(const complex<double>& x);

    double bsh21m1_outside_2_coeff(int i, int j);
    double bsh21m1_outside_1_coeff(int i);


    complex<double> mysign(const complex<double>& x);
    
    
    int find_region(const complex<double>& x);
    
    complex<double> cli2(const complex<double>& z);
    complex<double> cli3(const complex<double>& z);
    complex<double> cli4(const complex<double>& z);
    complex<double> cli4_sbc(const complex<double>& z);
    complex<double> cli4_sbc_2(const complex<double>& z);
    complex<double> ch2m2(const complex<double>& z);
    complex<double> ch21m1(const complex<double>& z);
    
    void check_index(const vector<int>& v);
    void check_index(int n1);
    int find_index_in_base_3(const vector<int>& vint);
}

#endif

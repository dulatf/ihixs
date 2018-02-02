#ifndef PSI_H
#define PSI_H

#include<complex>
#include<vector>
using namespace std;


complex<double> psi(const complex<double>& z);

class PsiFinder
{
private:
    vector<double> imag_z;
    vector<double> real_psi;
    vector<double> imag_psi;
    vector<double> real_psi_p_1;
    vector<double> imag_psi_p_1;
    vector<double> real_psi_p_2;
    vector<double> imag_psi_p_2;
    
    static PsiFinder *s_instance;
    PsiFinder(){};// private constructor
public:
    complex<double> get(const complex<double>&);
    

    static PsiFinder *psi()
    {
        if (!s_instance)
            s_instance = new PsiFinder;
        return s_instance;
    }
};

#endif

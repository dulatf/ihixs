#ifndef LUMINOSITY_H
#define LUMINOSITY_H

/**
 *
 * \file    luminosity.h
 * \ingroup tools
 * \author  Achilleas Lazopoulos
 * \date    February 2015
 *
 */

#include <string>
using namespace std;

#include "LHAPDF/LHAPDF.h"


/**
 *
 * \class   Luminosity
 * \ingroup tools
 * \brief   Interface with LHAPDFs. It initializes a PDF member (one of the many of a specific PDF grid). Then the user can add luminosity pairs (initial states) with a potential constant factor that depends on the pair. Once configured, the class provides L(x1,x2,muf) = Sum_i {x1*f1_i(x1,muf) * x2*f2_i(x2,muf) * c_i}
 * \todo
 *
 */

class Luminosity{
private:
    /// \name Private data members
    /// @{
    /// \brief vector of pointers to a whole PDF member (includes all flavors)
    vector<LHAPDF::PDF*> _pdf;///<

    /// \brief each pair is a different initial state
    vector<pair<int,int> > _pairs;

    /// \brief potential coefficient that is different for every initial state pair
    vector<double> _coeff;

    /// \brief technical cutoff: we don't allow bjorken xs to go closer to the edges than this cutoff.
    static const double _almost_zero;
    /// @}
public:
    /// \name Constructors and destructors
    /// @{
    Luminosity(const string& gridname);
    Luminosity(const string& gridname,int imember);

    ~Luminosity(){
            for(int i = 0; i < _pdf.size(); ++i)
                 delete _pdf[i];
            }
    ///< we are responsible for deleting the LHAPDF pointer

    ///@}

    /// \name Input/Output
    /// @{

    /// \brief adds a pair of initial states by number id
    void addPair(int left,int right);
    /// \brief adds a pair of initial states by number id, *and* a constant. The constant multiplies the luminosity contribution from this pair of initial states. (Constant coeff defaults to 1.0)

    void addPair(int left,int right,const double& c);
    /// \brief clears the class data
    /// \todo where is this used? Is it really useful?

    void clear_pairs(){_pairs.clear();_coeff.clear();}
    /// \brief returns the full luminosity Sum_i {x1*f1_i(x1,muf) * x2*f2_i(x2,muf) * c_i}

    double give(const double& x1,const double& x2,const double& muf);
    /// \brief returns the luminosity with pdfs set to _pdf[0] which might or might not be the central pdf member of the set. If the Luminosity was initialized through the Luminosity(const string& gridname,int imember), then the _pdf[0] points to the imember, otherwise it points to whatever member is the 0th member of this set according to LHAPDF classification (most of the times this is the central set, if there is one).

    vector<double> luminosity_vector(const double& x1,const double& x2,const double& muf);
    /// \brief returns a vector of lumis for each pdf member in the _pdf vector
    /// \warning it makes little sense to use this with the Luminosity(const string& gridname,int imember) constructor
    
    double give_specific(const double& x1,const double& x2,const double& muf, int imember);
    /// \brief lumi for the imember in the _pdf vector
    /// \warning it makes little sense to use this with the Luminosity(const string& gridname,int imember) constructor
    
    
    double alpha_s_at_mz(){return _pdf[0]->alphasQ(91.1876);}
    /// \brief alpha_s_at_mz
    /// \warning mz is at its nominal pdg value here - might not be in sync with model
    double as_at(const double mur);

    /// \brief sample and return pdf xf
    double give_single(const double &x, const double & muf, int flavor);


    /// \brief return the size of _pdf
    int Size(){return _pdf.size();}
    ///@}

};




#endif

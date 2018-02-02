#include <strstream>
#include <fstream>
#include "gauss_grid.h"
#include "cuba.h"
using namespace std;

/////////// legendre / contour grid ///////////
namespace Legendre{
  const unsigned int Order = 20;
  double Positions[Order] = {-0.99312859918509492478612238847132027822264713090166,-0.96397192727791379126766613119727722191206032780619,-0.91223442825132590586775244120329811304918479742369,-0.83911697182221882339452906170152068532962936506564,-0.74633190646015079261430507035564159031073067956918,-0.63605368072651502545283669622628593674338911679937,-0.51086700195082709800436405095525099842549132920243,-0.37370608871541956067254817702492723739574632170568,-0.22778585114164507808049619536857462474308893768293,-0.076526521133497333754640409398838211004796266813498,0.076526521133497333754640409398838211004796266813498,0.22778585114164507808049619536857462474308893768293,0.37370608871541956067254817702492723739574632170568,0.51086700195082709800436405095525099842549132920243,0.63605368072651502545283669622628593674338911679937,0.74633190646015079261430507035564159031073067956918,0.83911697182221882339452906170152068532962936506564,0.91223442825132590586775244120329811304918479742369,0.96397192727791379126766613119727722191206032780619,0.99312859918509492478612238847132027822264713090166};
  double Weights[Order] = {0.0176140071391521183118619623518528163621431,0.0406014298003869413310399522749321098790906,0.0626720483341090635695065351870416063516011,0.0832767415767047487247581432220462061001778,0.10193011981724043503675013548034987616669166,0.11819453196151841731237737771138228700504122,0.131688638449176626898494499748163134916110511,0.1420961093183820513292983250671649330345154134,0.149172986472603746787828737001969436692679904081,0.1527533871307258506980843319550975934919486451124,0.1527533871307258506980843319550975934919486451124,0.149172986472603746787828737001969436692679904081,0.1420961093183820513292983250671649330345154134,0.131688638449176626898494499748163134916110511,0.11819453196151841731237737771138228700504122,0.10193011981724043503675013548034987616669166,0.0832767415767047487247581432220462061001778,0.0626720483341090635695065351870416063516011,0.0406014298003869413310399522749321098790906,0.0176140071391521183118619623518528163621431};
}
namespace Contour{
  const unsigned int Order = 37;
  double Subdivisions[Order] = {0.,0.5,1.,1.5,2.,2.5,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125};
  double CPos = 2.5;
}
///////// end legendre / contour grid /////////


GaussGrid::GaussGrid(const UserInterface& UI, const string &gridPath) : _gridPath(gridPath){
    _verbose = UI.giveBool("with_resummation_info");
    _numIntegralSamples = 2000000;
    _lumi = new Luminosity(UI.giveString("pdf_set"),UI.giveInt("pdf_member"));
    _scale = UI.giveDouble("muf");
    _numContourPoints = Legendre::Order * (Contour::Order-1);
    _subdivisions.resize(_numContourPoints);
    _tau = pow(UI.giveDouble("m_higgs")/UI.giveInt("Etot"),2);
    if (_verbose)
        cout << "[GaussGrid] Generating subdivisions..." << endl;
    generate_subdivisions();
    if (_verbose)
        cout << "[GaussGrid] Subdivisions done..." << endl;
    bool gridExists = check_grid(gridPath);
    if(gridExists){
        cout << "[GaussGrid] " << "Grid " << gridPath << " exists." << endl;
        load_grid();
        export_grid("test.txt");
    }
    else{
        cout << "[GaussGrid] " << "Grid " << gridPath << " does not exixt." << endl;
        cout << "[GaussGrid] " << "We generate a new one. This can take some time, please be patient" << endl;
        generate_grid(UI);
        export_grid("test.txt");
    }
}
GaussGrid::~GaussGrid(){
  delete _lumi;
}

double GaussGrid::DoIntegral(const function<complex<double> (complex<double>)> &integrand){
  complex<double> result = 0.0;
  vector<complex<double> >::iterator itn;
  vector<pair<double, double> >::iterator itp;
  cout << setprecision(16);
  for(itn = _moments.begin(), itp = _subdivisions.begin(); itn != _moments.end() && itp != _subdivisions.end(); itn++, itp++){
    complex<double> cn(Contour::CPos,itp->first);
     // cout<<cn<<","<<endl;
    complex<double> nt = integrand(cn);
    complex<double> cc = itp->second * pow(*itn,2) * nt;
    result += cc;
  }
  return result.real()/consts::Pi;
}


bool GaussGrid::check_grid(const string &gridPath){
  ifstream f(gridPath);
  return f.good();
}

void GaussGrid::load_grid(void){
    if (_verbose)
        cout << "[GaussGrid] " << "Loading grid from " << _gridPath << "" << endl;
  ifstream fi(_gridPath,ios::binary);
  double tmp[2];
  _moments.resize(_numContourPoints);
  for(unsigned int i=0;i < _numContourPoints;i++){
    fi.read((char*)tmp,sizeof(double)*2);
    _moments[i] = complex<double>(tmp[0],tmp[1]);
  }
}

void GaussGrid::export_grid(const string &path){
  ofstream of(path);
  of << "{" << fixed;
  vector<complex<double> >::iterator it;
  vector<pair<double, double> >::iterator pit;
  //for(vector<complex<double> >::iterator it=_moments.begin();it!=_moments.end();it++){
  for(it = begin(_moments), pit=begin(_subdivisions);it!=end(_moments)&&pit!=end(_subdivisions);it++,pit++){
  //for(auto it : _moments, auto pit : _subdivisions){
    complex<double> cn(Contour::CPos,pit->first);
    of << real(cn) << " + I * " << imag(cn) << ",\t" << real(*it) << " + I * " << imag(*it) << ",\t" << pit->second;
     if(it+1 != _moments.end())
       of << ",\n";
  }
  of << "}";
}
void GaussGrid::store_grid(void){
  ofstream of(_gridPath,ios::binary|ios::out|ios::trunc);
  for(vector<complex<double> >::iterator it=_moments.begin();it!=_moments.end();it++){
    double re = real(*it);
    double im = imag(*it);
    of.write((char*)&re,sizeof(double));
    of.write((char*)&im,sizeof(double));
  }
}
void GaussGrid::generate_grid(const UserInterface& UI){
    if (_verbose){
        cout << "[GaussGrid] " << "Generating a new grid" << endl;
        cout << "[GaussGrid] " << "Pdf used: " << UI.giveString("pdf_set") << "/" << UI.giveInt("pdf_member") << endl;
        cout << "[GaussGrid] " << "Factorization scale: " << UI.giveDouble("muf") << endl;
        cout << "[GaussGrid] Computing the grid..." << endl;
    }
    compute_moments();
    if (_verbose)
        cout << "[GaussGrid] Finished computing the grid" << endl;
    store_grid();
}

complex<double> GaussGrid::mellin_integral(const complex<double> &n){
  double dx = 1.0 / (double)_numIntegralSamples;
  complex<double> acc(0.0,0.0);
  for(unsigned int i = 0; i < _numIntegralSamples;i++){
    double x = i * dx;
    double pdf = _lumi->give_single(x,_scale,QCD::g);
    acc += pow(x,n-2.0) * pdf;
  }
  return acc * dx;
}

/*void GaussGrid::compute_moments(void){
  _moments.resize(_numContourPoints);
   #pragma omp parallel for
  for(unsigned int i = 0; i < _numContourPoints;i++){
    cout << ".";flush(cout);
    if(i%5==0){
      cout << "\r\x1b[K";
      flush(cout);
    }
    complex<double>cn(Contour::CPos,_subdivisions[i].first);
    _moments[i] = mellin_integral(cn);
  }
  cout << "\n";
}*/
/*struct GIntParam{
    complex<double> cn;
    Luminosity *lumi;
    double scale;
};
int integrand(const int *ndim, const double *x, const int *ncomp, double *f, void *userdata){
    GIntParam *gp = (GIntParam*)userdata;
    double pdf = gp->lumi->give_single(x[0],gp->scale,QCD::g);
    complex<double> rr = pow(x[0],gp->cn-2.0)*pdf;
    f[0] = real(rr);
    f[1] = imag(rr);
    return 0;
}
void GaussGrid::compute_moments(void){
    _moments.resize(_numContourPoints);
    int ndim=1;
    int ncomp=2;
    double epsrel = 1e-6;
    double epsabs = 1e-12;
    int flags=0;
    int seed=0;
    int mineval=1000;
    int maxeval=1000000;
    int nstart=100;
    int nstep=100;
    int nbatch=5000;
    
    for(unsigned int i=0;i<_numContourPoints;i++){
        cout << '.';flush(cout);
        if(i%10==0){
            cout << "\r\x1b[K";
            flush(cout);
        }
        GIntParam gp{{Contour::CPos,_subdivisions[i].first},_lumi,_scale};
        double val[2], error[2], prob[2];
        int neval, fail;
        Vegas(ndim,ncomp,integrand,&gp,epsrel,epsabs,flags,seed,mineval,maxeval,nstart,nstep,nbatch,0,nullptr,&neval,&fail,val,error,prob);
        _moments[i] = complex<double>(val[0],val[1]);
    }
    cout << endl;
}
*/
/*
 struct GIntParam{
    complex<double> *cns;
    Luminosity *lumi;
    double scale;
};
int integrand(const int *ndim, const double *x, const int *ncomp, double *f, void *userdata){
    GIntParam *gp = (GIntParam*)userdata;
    double pdf = gp->lumi->give_single(x[0],gp->scale,QCD::g);
    for(int i = 0; i < *ncomp;i+=2){;
        complex<double> rr = pow(x[0],gp->cns[i/2]-2.0)*pdf;

        f[i] = real(rr);
        f[i+1] = imag(rr);
    }
    return 0;
}
void GaussGrid::compute_moments(void){
    _moments.resize(_numContourPoints);
    cout << "We're gonna try this with " << 2*_numContourPoints <<  " components now."  << endl;
    int ndim=1;
    int ncomp=2*_numContourPoints;
    double epsrel = 1e-3;
    double epsabs = 0;///1e-7;
    int flags=0;
    int seed=0;
    int mineval=1000;
    int maxeval=10000000;
    int nstart=1000;
    int nstep=1000;
    int nbatch=5000;
    
    complex<double> *cns = new complex<double>[_numContourPoints];
    for(unsigned int i=0;i<_numContourPoints;i++){
        cns[i] = complex<double>{Contour::CPos,_subdivisions[i].first};
    }
    GIntParam gp{cns,_lumi,_scale};
    double *val = new double[ncomp];
    double *error = new double[ncomp];
    double *prob = new double[ncomp];
    int neval, fail;
    Vegas(ndim,ncomp,integrand,&gp,epsrel,epsabs,flags,seed,mineval,maxeval,nstart,nstep,nbatch,0,nullptr,&neval,&fail,val,error,prob);
     
    cout << maxeval << " vs "  << neval <<  " -> done" << endl;

    for(unsigned int i=0;i<_numContourPoints;i++){
            _moments[i] = complex<double>{val[2*i+0],val[2*i+1]};
    }
    delete [] val;
    delete [] error;
    delete [] prob;
    delete [] cns;
}

*/
struct GIntParam{
    complex<double> *cns;
    Luminosity *lumi;
    double scale;
    int start;
    double tau;
};

int integrand(const int *ndim, const double *x, const int *ncomp, double *f, void *userdata){
    GIntParam *gp = (GIntParam*)userdata;
    double xx = x[0];//gp->tau * (1-gp->tau)*x[0];
    double pdf = gp->lumi->give_single(xx,gp->scale,QCD::g);
    for(int i = gp->start; i < (gp->start+*ncomp);i+=2){
        complex<double> rr = pow(xx,gp->cns[i/2]-2.0)*pdf;

        f[i-gp->start] = real(rr);
        f[i+1-gp->start] = imag(rr);
    }
    return 0;
}
void GaussGrid::compute_moments(void){
    _moments.resize(_numContourPoints);
    if (_verbose) cout << "We're gonna try this with " << 2*_numContourPoints <<  " components now."  << endl;
    int blockSize = 5;
    int numBlocks = (_numContourPoints / blockSize);
    int ndim=1;
    int ncomp=2*blockSize;
    double epsrel = 1e-3;
    double epsabs = 0;///1e-7;
    int flags=0;
    int seed=0;
    int mineval=1000;
    int maxeval=5000000;
    int nstart=10000;
    int nstep=100;
    int nbatch=500;
    int spin[1];
    spin[0]= -1;// see cuba doc for the meaning of this

    
    complex<double> *cns = new complex<double>[_numContourPoints];
    for(unsigned int i=0;i<_numContourPoints;i++){
        cns[i] = complex<double>{Contour::CPos,_subdivisions[i].first};
    }
    GIntParam gp{cns,_lumi,_scale,0,_tau};
    double *val = new double[ncomp];
    double *error = new double[ncomp];
    double *prob = new double[ncomp];
    for(int i=0;i<numBlocks;i++){
        int neval, fail;
        gp.start = blockSize * i * 2;
        if (_verbose) cout << "Iteration #" << i << " start is at " << gp.start << endl;
        Vegas(ndim,ncomp,integrand,&gp,1,epsrel,epsabs,flags,seed,mineval,maxeval,nstart,nstep,nbatch,0,nullptr,spin,&neval,&fail,val,error,prob);
        
     
        if (_verbose) cout << maxeval << " vs "  << neval <<  " -> done" << endl;
        for(unsigned int j=0;j<ncomp;j++){
            if (_verbose) {
            cout << ">>\t" << val[j] << " +/- " << error[j]
                << " : " << prob[j]
                << " \t error/val x 100% = "<<error[j]/val[j]*100.<<endl;
            }
        }
        for(unsigned int j=0;j<blockSize;j++){
            _moments[blockSize*i+j] = complex<double>{val[2*j+0],val[2*j+1]};
        }
        if (_verbose) cout << "end of Iteration #" << i << endl;
        else cout << ".";
    }

    delete [] val;
    delete [] error;
    delete [] prob;
    delete [] cns;
}



void GaussGrid::generate_subdivisions(void){
  for(unsigned int sd = 0; sd < Contour::Order-1; sd++){
    double a = Contour::Subdivisions[sd];
    double b = Contour::Subdivisions[sd+1];
    auto hx = [a,b](double y) -> double { return a + (b-a) * (1.0 + y) * 0.5; };
    double jac = 0.5 * (b-a);
      if (_verbose)
          cout << "[GaussGrid] Sampling interval (" << a << ", " << b << ")\r" ;
    for(unsigned int l = 0; l < Legendre::Order; l++){
      _subdivisions[sd*Legendre::Order + l] = make_pair(hx(Legendre::Positions[l]),jac * Legendre::Weights[l]);
    }
  }
  cout << endl;
}

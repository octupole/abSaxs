/*
 * Funktionell.h
 *
 *  Created on: Dec 13, 2018
 *      Author: marchi
 */

#ifndef ABINITIO_FUNKTIONELL_H_
#define ABINITIO_FUNKTIONELL_H_
#include <map>
#include <functional>
#include <utility>
#include <complex>
#include "MyUtilClass.h"
#include "Ftypedefs.h"
#include "RhoSaxs.h"
#include "SaxsData.h"
using namespace MATRIX;
using namespace DVECT;

#include "Array.h"
using Complex=std::complex<double>;

using uint=unsigned int;
using std::stringstream;
using std::map;
using std::tuple;
using namespace Array;
using namespace abinit;
namespace Funkll {
class Funktionell {
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	uint nx{0},ny{0},nz{0},nzp{0};
	uint Nx{0},Ny{0},Nz{0};
	Matrix CO,OC,co,oc;
	double dq{0.002};
	double qcut{1.0};
	map<size_t,std::pair<double,size_t>> Iq_c;
	map<size_t,double> Iq_exp;
	map<size_t,vector<vector<int>>> mapIdx;
	map<size_t,vector<DDvect<double>>> mapVec;

	double Scaling{0};
	double Gscaling{0};
	void setUpFirst(SaxsData *);
	array3<Complex> Modulus(array3<Complex> &);
	double Par{1.0e8};
public:
	Funktionell()=delete;
	Funktionell(RhoSaxs *,RhoSaxs *,SaxsData *);
	double Energy(array3<Complex> &);
	double GradPar(){return Gscaling;}
	double Deviate2(array3<Complex> &);
	virtual ~Funktionell();
};

} /* namespace Funkll */

#endif /* ABINITIO_FUNKTIONELL_H_ */
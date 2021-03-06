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

#include "abRhoSaxs.h"
#include "MyUtilClass.h"
#include "Ftypedefs.h"
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

namespace Funkll {
class Funktionell {
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	uint nx{0},ny{0},nz{0},nzp{0};
	uint Nx{0},Ny{0},Nz{0};
	Matrix CO,OC,co,oc;
	double dq{0.002};
	double qcut{1.0};
	double qmin{0};
	double Rd{0},Pot{1.0e6};
	double Rge{0};
	map<size_t,double> Iq_c;
	map<size_t,double> Iq_exp;
	map<size_t,vector<vector<int>>> mapIdx;
	map<size_t,vector<DDvect<double>>> mapVec;

	void setUpFirst(SaxsData *);
	array3<Complex> Modulus(array3<Complex> &);
	double Par{100.0};
	double A_1{1},A_0{0};
public:
	Funktionell()=delete;
	Funktionell(abInitioRho::RhoSaxs *,abInitioRho::RhoSaxs *,SaxsData *);
	double EnergyQ(double &,double &,array3<Complex> &);
	double EnergyQ(double &,double &,array3<Complex> &, array3<Complex> &);

	void ComputeIqc(array3<Complex> &);
	double EnergyR(double &,array3<double> &, array3<double> &);
	map<size_t,double> & getIqc(){return Iq_c;};
	map<size_t,double> & getIqe(){return Iq_exp;};
	map<size_t,vector<vector<int>>> & getIdx(){return mapIdx;};
	void scaleIq(array3<Complex> &);
	void Write();
	double getMydq(){return dq;};
	virtual ~Funktionell();
};

} /* namespace Funkll */

#endif /* ABINITIO_FUNKTIONELL_H_ */

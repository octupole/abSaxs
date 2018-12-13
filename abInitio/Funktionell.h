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
#include <tuple>
#include <complex>
#include "MyUtilClass.h"
#include "Ftypedefs.h"
#include "RhoSaxs.h"
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
	map<int,tuple<double,double>> Iq_c;
	array3<Complex> Grad;
public:
	Funktionell()=delete;
	Funktionell(RhoSaxs &,RhoSaxs &);
	double Deviate(array3<Complex> &);
	virtual ~Funktionell();
};

} /* namespace Funkll */

#endif /* ABINITIO_FUNKTIONELL_H_ */

/*
 * ABSaxs.h
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#ifndef ABINITIO_ABSAXS_H_
#define ABINITIO_ABSAXS_H_
#include "SaxsData.h"
#include "MyUtilClass.h"
#include "Ftypedefs.h"
#include "Prcfft3d.h"
#include "Pcrfft3d.h"
#include "Funktionell.h"
#include <mpi.h>

#include <sstream>
#include <map>
#include <functional>
#include <tuple>

#include "abRhoSaxs.h"
using namespace MATRIX;
using namespace DVECT;
using uint=unsigned int;
using std::stringstream;
using std::map;
using std::tuple;
namespace abinit{
const double CellParam{3.5};

class ABSaxs {
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	SaxsData * Exp{nullptr};
	SaxsData * Calc{nullptr};
	Matrix CO{0},co{0};
	Matrix OC{0},oc{0};
	vector<uint> grid_a{3};
	vector<uint> grid_b{3};
	uint nx{0},ny{0},nz{0};
	uint Nx{0},Ny{0},Nz{0};
	uint nzp{0};
	double dx{0},dy{0},dz{0},dvol{0};
	double unitsR{10.0};
	double unitsQ{1.0/unitsR};

	double dq{0.002};
	double qcut{1.0};
	double SuperCell{1.0};
	std::function<void(double)> metric=[this](double R,double c=CellParam)
			{co[XX][XX]=R*c;co[YY][YY]=R*c;co[ZZ][ZZ]=R*c;CO=SuperCell*co;};
	abInitioRho::RhoSaxs * Rho_in{nullptr}, * Rho_s{nullptr};
	Pfftwpp::Prcfft3d * Forward3{nullptr};
	Pfftwpp::Pcrfft3d * Backward3{nullptr};
	array3<Complex> F_k;
	array3<double> F_r;
	array3<double> GradR;

	array3<Complex> I_k;
	array3<Complex> Modulus(array3<Complex> &);
	Funkll::Funktionell * myFuncx{nullptr};
	void __qhistogram();
	double scalePlot{0};
	double Gscale{0};

public:
	ABSaxs()=delete;
	ABSaxs(uint, uint, uint, double=1.0);
	void getMetrics(Matrix & bCO,Matrix & bco){bCO=CO;bco=co;}
	void setUpCell(SaxsData *);
	void setUpRho(SaxsData *);
	void setUpRho(SaxsData *,array3<double> &);

	void Run();
	virtual void Minimize();
	void testGradient();
	void testMinim();
	void minLBFGS(const real_1d_array &, double &, real_1d_array &, void *);
	virtual ~ABSaxs();
};
}
#endif /* ABINITIO_ABSAXS_H_ */

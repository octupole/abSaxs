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
#include "RhoSaxs.h"
#include <functional>
using namespace MATRIX;
using namespace DVECT;
using uint=unsigned int;

namespace abinit{
const double CellParam{3.5};
class ABSaxs {
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	SaxsData * Exp{nullptr};
	SaxsData * Calc{nullptr};
	Matrix CO{0},co{0};
	vector<uint> grid_a{3};
	vector<uint> grid_b{3};
	double SuperCell{1.0};
	std::function<void(double)> metric=[this](double R,double c=CellParam)
			{co[XX][XX]=R*c;co[YY][YY]=R*c;co[ZZ][ZZ]=R*c;CO=SuperCell*co;};
	RhoSaxs Rho_in,Rho_s;
public:
	ABSaxs()=delete;
	ABSaxs(uint, uint, uint, double);
	void setUp(SaxsData *);
	void setUp();
	void Run();
	virtual ~ABSaxs();
};
}
#endif /* ABINITIO_ABSAXS_H_ */

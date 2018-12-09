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

#include "Grid.h"
#include <functional>
using namespace MATRIX;
using namespace DVECT;


namespace abinit{
class ABSaxs {
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	SaxsData * Exp{nullptr};
	SaxsData * Calc{nullptr};
	Matrix CO;
	std::function<void(double)> metric=[this](double R,double c=6.0)
			{CO[XX][XX]=R*c;CO[YY][YY]=R*c;CO[ZZ][ZZ]=R*c;};
public:
	ABSaxs();
	void setUp(SaxsData *);
	void setUp();
	void Run();
	virtual ~ABSaxs();
};
}
#endif /* ABINITIO_ABSAXS_H_ */

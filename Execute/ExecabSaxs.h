/*
 * ExecabSaxs.h
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */

#ifndef EXECUTE_EXECABSAXS_H_
#define EXECUTE_EXECABSAXS_H_

#include "TrjRead.h"
#include "MyUtilClass.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include "SaxsData.h"
#include "ABSaxs.h"
#include "Saxs.h"
#include "abRhoSaxs.h"
#include "RhoSaxs.h"
#include "RhoSaxsLI.h"
#include "Array.h"
#include "Atoms.h"
#include "Metric.h"
#include "Topol.h"

using std::vector;
using Matrix=MATRIX::MMatrix<double>;
using MAtoms=Atoms<double>;
class ExecabSaxs {
	Matrix CO,co;
	size_t nx{1},ny{1},nz{1},Nx{1},Ny{1},Nz{1};
	ifstream * ExIn;
	ofstream * ExOut;
	SaxsData * Exp{nullptr};
	DensMode myDens;
	int myDensAvg{0};

	abinit::ABSaxs * myABSaxs{nullptr};
	void expSaxs();
	Saxs * MySaxs{nullptr}; ///< Pointer to the Saxs class to be instantiated
	RhoSaxs * Rho_ex{nullptr};
	Topol_NS::Topol * Top{nullptr};
	void __Saxs(MAtoms *);
	void __SuperCell();
public:
	ExecabSaxs(trj::TrjRead &,Topol_NS::Topol & );
	ExecabSaxs()=delete;
	void Run_abSaxs(MAtoms *);
	virtual ~ExecabSaxs();
};

#endif /* EXECUTE_EXECABSAXS_H_ */

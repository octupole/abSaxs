/*
 * Spline1DInterpolant.h
 *
 *  Created on: Mar 22, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_SPLINE1DINTERPOLANT_H_
#define TRJLIB_SPLINE1DINTERPOLANT_H_
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <set>
#include "stdafx.h"
#include "interpolation.h"
#include <map>
#include "histograms.hpp"

using namespace alglib;
using namespace std::placeholders;
using std::bind;

using std::function;
using std::vector;
using std::set;
using std::map;
namespace Spline1D {

class Spline1DInterpolant {
	real_1d_array x,y;
	spline1dinterpolant s;
	double Dq{0},myUnits{1.0},cutoff{0};

public:
	Spline1DInterpolant()=delete;
	Spline1DInterpolant(vector<double>,vector<double>);
	Spline1DInterpolant(double,map<size_t,double> &);
	Spline1DInterpolant(const vector<vector<double> >& );
	Spline1DInterpolant(Histogram1D *, double,double);
	double lowLimit()const {return x[0];}
	double operator()(double);
	Spline1DInterpolant & operator-=(const Spline1DInterpolant &);
	function<double(const double)> Spline;
	virtual ~Spline1DInterpolant();
	friend ostream & operator<<(ofstream &,Spline1DInterpolant &);
};

} /* namespace Spline1D */

#endif /* TRJLIB_SPLINE1DINTERPOLANT_H_ */

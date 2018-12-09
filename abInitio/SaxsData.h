/*
 * SaxsData.h
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#ifndef UTILS_SAXSDATA_H_
#define UTILS_SAXSDATA_H_
#include <vector>
#include <tuple>
#include <fstream>
#include <iostream>
#include "interpolation.h"
#include "SGSmooth.h"
#include <functional>

using std::vector;
using std::tuple;
using std::ostream;
using std::ifstream;
using std::endl;
using std::cout;
using namespace alglib;

class SaxsData {
	vector<tuple<double,double>> x;
	double Rg{0};
	void Generate(vector<double> & ,vector<double> &);
public:
	std::function<double(void)>  Rd{[this](){return sqrt((5.0/3.0)*Rg*Rg);}};
	double getRg(){return Rg;}
	SaxsData(vector<double> &,vector<double> &);
	SaxsData()=delete;
	tuple<double,double> & operator[](size_t);
	virtual ~SaxsData();
	friend ostream & operator<<(ostream &, SaxsData &);

};

#endif /* UTILS_SAXSDATA_H_ */

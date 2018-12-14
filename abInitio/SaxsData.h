/*
 * SaxsData.h
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#ifndef UTILS_SAXSDATA_H_
#define UTILS_SAXSDATA_H_
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include "interpolation.h"
#include "SGSmooth.h"
#include <functional>
#include <algorithm>

using std::vector;
using std::ostream;
using std::ifstream;
using std::endl;
using std::cout;
using namespace alglib;

class SaxsData {
	vector<std::pair<double,double>> x;
	double Rg{0};
	void Generate(vector<double> & ,vector<double> &);
public:
	std::function<double(void)>  Rd{[this](){return sqrt((5.0/3.0)*Rg*Rg);}};
	double getRg(){return Rg;}
	SaxsData(vector<double> &,vector<double> &);
	SaxsData()=delete;
	std::pair<double,double> & operator[](size_t);
	vector<std::pair<double,double>> & gIq_exp(){return x;}
	virtual ~SaxsData();
	friend ostream & operator<<(ostream &, SaxsData &);

};

#endif /* UTILS_SAXSDATA_H_ */

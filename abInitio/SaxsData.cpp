/*
 * SaxsData.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#include "SaxsData.h"

SaxsData::SaxsData(vector<double> & x,vector<double> & y) {
//	this->x=vector<tuple<double,double>>(x.size());
//	for(size_t o{0};o<x.size();o++){
//		this->x[o]=std::make_tuple(x[o],y[o]);
//	}
	this->Generate(x,y);
}

void SaxsData::Generate(vector<double> & x,vector<double> & y){
	const int w{100},deg{6};
	double qcut{0.06};
	real_1d_array xd,yd,a2;
	ae_int_t basisFuncs{2};
	double t{2},v;
	ae_int_t info;
	barycentricinterpolant p;
	polynomialfitreport rep;
	size_t regrMax{0};

	vector<double> y_s;
	y_s=sg_smooth(y,w,deg);
	this->x=vector<tuple<double,double>>(x.size());
	for(size_t o{0};o<x.size();o++){
		this->x[o]=std::make_tuple(x[o],y_s[o]);
	}

	for( auto x0: x){
		if(x0 <= qcut) regrMax++;
	}
	regrMax-=2;
	size_t mm{regrMax};

	xd.setlength(mm);
	yd.setlength(mm);
	for(int o=0;o<mm;o++){
		xd[o]=x[o+2]*x[o+2];
		yd[o]=log(y_s[o+2]);
	}
	polynomialfit(xd, yd, basisFuncs, info, p, rep);

	v = barycentriccalc(p, t);

	polynomialbar2pow(p, a2);
	for(auto o{0};o<a2.length();o++)
		std::cout << a2[o]<< " ";
	Rg=sqrt(-3.0*a2[1]);
	cout <<" Rg "<< Rg<< " " << Rd() <<endl;

//	cout << *this;
}

tuple<double,double> & SaxsData::operator [](size_t N){
	return x.at(N);
}

SaxsData::~SaxsData() {
	// TODO Auto-generated destructor stub
}

ostream & operator<<(std::ostream & fout, SaxsData & y){
	fout << "# Gyration ratio Rg = "<< y.Rg << std::endl;
	for(auto tp: y.x)
		fout << std::get<0>(tp)<<" " << std::get<1>(tp)<< std::endl;

	return fout;
}


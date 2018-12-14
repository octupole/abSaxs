/*
 * SaxsData.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#include "SaxsData.h"

SaxsData::SaxsData(vector<double> & x,vector<double> & y) {

	this->Generate(x,y);
}

void SaxsData::Generate(vector<double> & x,vector<double> & y){
	const int W{100},DEG{6};
	const int SEGM{10},SKIP{5};
	double qcut{0.07};
	ae_int_t basisFuncs{2};
	ae_int_t info;
	barycentricinterpolant p;
	polynomialfitreport rep;
	size_t regrMax{0};

	vector<double> y_s;
	y_s=sg_smooth(y,W,DEG);

	this->x=vector<std::pair<double,double>>(x.size());
	for(size_t o{0};o<x.size();o++){
		this->x[o].first=x[o];
		this->x[o].second=y_s[o];
	}
	for( auto x0: x){
		if(x0 <= qcut) regrMax++;
	}
	size_t mm{regrMax};

	int M{(int) mm-SEGM};
	int Mtime{M/SKIP};

	vector<double> Rg_s,err0;
	for(int o{0};o< M;o+=SKIP){
		double Rg{0};
		ae_int_t basisFuncs{2};
		double t{2},v;
		ae_int_t info;
		barycentricinterpolant p;
		polynomialfitreport rep;
		real_1d_array a2;
		real_1d_array xd,yd;
		int iBeg{o},iEnd{o+SEGM};
		int n{0};
		xd.setlength(iEnd-iBeg+1);
		yd.setlength(iEnd-iBeg+1);
		for(int p{iBeg};p<=iEnd;p++,n++){
			xd[n]=x[p]*x[p];
			yd[n]=log(y_s[p]);
		}
		polynomialfit(xd, yd, basisFuncs, info, p, rep);

		v = barycentriccalc(p, t);

		polynomialbar2pow(p, a2);

		Rg=sqrt(-3.0*a2[1]);
		Rg_s.push_back(Rg);
		err0.push_back(rep.avgerror);
	}
	auto it=std::min_element(err0.begin(),err0.end());
	Rg=Rg_s[std::distance(err0.begin(),it)];

}

std::pair<double,double> & SaxsData::operator [](size_t N){
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


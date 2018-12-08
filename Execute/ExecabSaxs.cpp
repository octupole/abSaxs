/*
 * ExecabSaxs.cpp
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */

#include "ExecabSaxs.h"

ExecabSaxs::ExecabSaxs(trj::TrjRead & MyIn) {
	ExIn=MyIn.gFin();
	ExOut=MyIn.gFout();
}
void ExecabSaxs::Run_abSaxs(){
	vector<double> x,y,y_s;
	double x0,y0,z0;
	while(*ExIn>>x0){
		*ExIn >> y0>>z0;
		x.push_back(x0);
		y.push_back(y0);
	}
	cout << y.size()<< endl;
	const int w{100},deg{6};
	y_s=sg_smooth(y,w,deg);
	for(size_t o{0};o<y_s.size();o++){
		*ExOut << x[o]<< " "<< y_s[o]<< " " << y[o] <<endl;
	}
}
ExecabSaxs::~ExecabSaxs() {
	// TODO Auto-generated destructor stub
}


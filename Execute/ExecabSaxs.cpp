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
	myABSaxs=new abinit::ABSaxs();
}
void ExecabSaxs::Run_abSaxs(){
	this->expSaxs();
	myABSaxs->setUp(Exp);
	myABSaxs->Run();
}
void ExecabSaxs::expSaxs(){
	vector<double> x,y,y_s;
	double x0,y0,z0;
	while(*ExIn>>x0){
		*ExIn >> y0>>z0;
		x.push_back(x0);
		y.push_back(y0);
	}
	Exp=new SaxsData(x,y);
}
ExecabSaxs::~ExecabSaxs() {
	// TODO Auto-generated destructor stub
}


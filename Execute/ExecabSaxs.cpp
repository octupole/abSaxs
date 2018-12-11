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
	myABSaxs=new abinit::ABSaxs(MyIn.gnnx(),MyIn.gnny(),MyIn.gnnz(),MyIn.gSupCell());
}
void ExecabSaxs::Run_abSaxs(){
	this->expSaxs();
	myABSaxs->setUp(Exp);
	myABSaxs->Run();
}
void ExecabSaxs::expSaxs(){
	vector<double> x,y,y_s;
	double x0,y0,z0;
	std::string line;
	try{
		while (getline(*ExIn, line)) {
			if(std::string{"#%"}.find(line.substr(0,1)) != std::string::npos)
				continue;
			std::istringstream iss{line};
			vector<string> field(std::istream_iterator<string>{iss},std::istream_iterator<string>());
			if(field.size()<2) throw string("At least two fields are needed in input data");
			std::stringstream{field[0]}>>x0;
			std::stringstream{field[1]}>>y0;
			x.push_back(x0);
			y.push_back(y0);
		}
	}catch(const string & s ){
		cout << s <<endl;
		exit(1);
	}
	Exp=new SaxsData(x,y);
}
ExecabSaxs::~ExecabSaxs() {
	// TODO Auto-generated destructor stub
}


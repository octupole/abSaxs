/*
 * ExecabSaxs.cpp
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */

#include "ExecabSaxs.h"

ExecabSaxs::ExecabSaxs(trj::TrjRead & MyIn, Topol_NS::Topol & Top0): Top{&Top0} {
	ExIn=MyIn.gFin();
	ExOut=MyIn.gFout();
	myABSaxs=new abinit::ABSaxs(MyIn.gnnx(),MyIn.gnny(),MyIn.gnnz(),MyIn.gSupCell());
	int MyOrder=MyIn.gMyOrder();
	double Myd=MyIn.gMyd();
	double Mycut=MyIn.gMyCut();
	if(MyOrder>1) Rho_ex=new RhoSaxsLI;
	else Rho_ex=new RhoSaxs;
	MySaxs=new Saxs(MyOrder,Myd,Mycut);
	nx=MyIn.gnnx();
	ny=MyIn.gnny();
	nz=MyIn.gnnz();
	MySaxs->Allocate(nx,ny,nz);
	auto fileDens=MyIn.gfiledens();
	myDens=MyIn.gModeCompute();
	myDensAvg=MyIn.gMyDensAvg();
	MySaxs->setDens(myDensAvg,fileDens,myDens);
}
void ExecabSaxs::__SuperCell(){
	MySaxs->setSuperCell0(CO[XX][XX]/co[XX][XX]);
	double SuperCell=CO[XX][XX];
	try{
		if(SuperCell < 0) throw string("pdb box does not have CRYST1 keyword set, cannot compute.");
	}catch(const string & s){cout << s <<endl;Finale::Finalize::Final();}
	MySaxs->setSuperCell(SuperCell);
	Nx=nx;
	Ny=ny;
	Nz=nz;
}

array3<double> & ExecabSaxs::__Saxs(MAtoms * atm){
	__SuperCell();
	MySaxs->Setup(atm->getIndx(),Top->get_atSFactor(),false,true);

	Rho_ex->Allocate(nx,ny,nz);
	MySaxs->ComputeDENS(Rho_ex,atm);
	ofstream fout_pdb("output.pdb",ios::out);
	atm->PrintAll(fout_pdb);
	fout_pdb.close();
	return MySaxs->getI_r();
}
void ExecabSaxs::Run_abSaxs(MAtoms * atm){
	atm->InertiaRot();
	this->expSaxs();
	if(atm){
		myABSaxs->setUpCell(Exp);
		myABSaxs->getMetrics(CO,co);
		Matrix COn=CO/10.0;
		Metric<double> Mt(COn);
		atm->setMT(Mt);
		array3<double> & F=__Saxs(atm);
		myABSaxs->setUpRho(Exp,F);
	} else{
		myABSaxs->setUpCell(Exp);
		myABSaxs->getMetrics(CO,co);
		myABSaxs->setUpRho(Exp);
	}
	//myABSaxs->testGradient();
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
	vector<std::pair<double,double>> & myExp=Exp->gIq_exp();
}
ExecabSaxs::~ExecabSaxs() {
	// TODO Auto-generated destructor stub
}


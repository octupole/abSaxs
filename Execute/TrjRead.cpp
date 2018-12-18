/*
 * TrjRead.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */

#include "TrjRead.h"
namespace trj {

TrjRead::TrjRead(int nv,char ** v): trjInput::trjInput(nv,v) {
	// TODO Auto-generated constructor stub

	string comm0(v[0]);
	size_t mypos=comm0.find_last_of("/")+1;
	size_t length=comm0.size()-mypos;
	string command=comm0.substr(mypos,length);
	string errmsg;
	string Usage="Usage:\t"+ command + "\n";
	vector<string> use=this->getUsage();
	vector<string> SelRes;
	string Reference;
	for(unsigned int n=0;n<use.size();n++)
		Usage+=use[n];
	Usage+="\n\t Default values in square brackets []\n";
	try{
		if(nv == 1) throw Usage;
		if(int m=this->bTest().size()) {
			errmsg=" Command(s) not found: ";
			for(int n=0;n<m;n++)
				errmsg+=this->bTest()[n]+"  ";
			errmsg+="\n"+Usage;
			throw errmsg;
			}
		}
	catch(const string & s){
		cout << s << endl;
		exit(1);
	}

}
void TrjRead::Input(){
	ifstream ftest;

	try{
		if(!inmap["-i"].empty()) {
			if(inmap["-i"].size() < 2) throw string("\n filename expected for " + inmap["-i"][0] + " option\n ");
			if(inmap["-i"].size() > 2) throw string("\n More than one entry for " + inmap["-i"][0] + " option \n");
			filein=inmap["-i"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();
			finx=new ifstream(filein.c_str(),ios::in);
		}
		if(!inmap["-o"].empty()) {
			if(inmap["-o"].size() < 2) throw string("\n filename expected for " + inmap["-o"][0] + " option \n");
			if(inmap["-o"].size() > 2) throw string("\n More than one entry for " + inmap["-o"][0] + " option \n");
			fileout=inmap["-o"][1];
			foutx=new ofstream(fileout.c_str(),ios::out);
		}
		if(!inmap["-supcell"].empty()) {
			int Nm=inmap["-supcell"].size();
			if(Nm == 2){
				stringstream(inmap["-supcell"][1])>> SuperCellSide;
			}
			else throw string(" One parameters is needed for  " + inmap["-supcell"][0] + " option ");
		}
		if(!inmap["-pdb"].empty()) {
			if(inmap["-pdb"].size() < 2) throw string("\n filename expected for " + inmap["-pdb"][0] + " option \n");
			if(inmap["-pdb"].size() > 2) throw string("\n More than one entry for " + inmap["-pdb"][0] + " option \n");
			filepdb=inmap["-pdb"][1];
			fpdb=new ifstream(filepdb.c_str(),ios::in);
			if(!fpdb) throw string("\n Cannot open " + filepdb + "!!\n");
		}
		if(!inmap["-dens"].empty()) {
			WhatToDo=Enums::ELDENS;
			int Nm=inmap["-dens"].size();
			std::stringstream ss;
			ss << Nm;
			ModeCompute(inmap["-dens"][1]);
			try{
				if(Nm == 3){
					stringstream(inmap["-dens"][2])>> MyOrder;
				} else if(Nm == 4){
					stringstream(inmap["-dens"][2])>> MyOrder;
					stringstream(inmap["-dens"][3])>> MyDensAvg;
				} else if(Nm > 4) throw string("\n Error: No. of parameters for "
						+ inmap["-grid"][0] + " exceeds three. It was "
						+ ss.str()+ " Abort!\n");
			}catch(const string & s){
				cout << s <<endl;
				exit(1);
			}
		}
		if(!inmap["-qhist"].empty()) {
			if(inmap["-qhist"].size() != 3) throw string(" Bin size and cutoff in q-space needed for  " + inmap["-qhist"][0] + " option ");
			stringstream(inmap["-qhist"][1])>> Myd;
			stringstream(inmap["-qhist"][2])>> Mycut;
		}

		if(!inmap["-grid"].empty()) {
			int Nm=inmap["-grid"].size();
			try{
				if(Nm == 2){
					stringstream(inmap["-grid"][1])>> nnx;
					nny=nnx; nnz=nnx;
				} else if(Nm == 4){
					stringstream(inmap["-grid"][1])>> nnx;
					stringstream(inmap["-grid"][2])>> nny;
					stringstream(inmap["-grid"][3])>> nnz;
				} else throw string("\n Warning: No Grid points provided " + inmap["-grid"][0] + " option. "
						+"Hope this is correct \n");
			}catch(const string & s){
				cout << s <<endl;
			}
		}

	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}

	gFin=finx;
	gFout=foutx;
	gFpdb=fpdb;
	if(WhatToDo==Enums::ELDENS){
#ifdef HAVE_VTK
		if(fileout.find_first_of(".") != std::string::npos){
			fileout=fileout.substr(0,fileout.find_first_of(".")+1)+"vts";
		} else{
			fileout=fileout+".vts";
		}
#else
		if(fileout.find_first_of(".") != std::string::npos){
			fileout=fileout.substr(0,fileout.find_first_of(".")+1)+"cvs";
		} else{
			fileout=fileout+".cvs";
		}
#endif

	}
}
TrjRead::~TrjRead() {

}

} /* namespace trj */

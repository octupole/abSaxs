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
		if(!inmap["-in"].empty()) {
			if(inmap["-in"].size() < 2) throw string("\n filename expected for " + inmap["-in"][0] + " option\n ");
			if(inmap["-in"].size() > 2) throw string("\n More than one entry for " + inmap["-in"][0] + " option \n");
			filein=inmap["-in"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();
			fin=ifstream(filein.c_str(),ios::in);
		}
		if(!inmap["-o"].empty()) {
			if(inmap["-o"].size() < 2) throw string("\n filename expected for " + inmap["-o"][0] + " option \n");
			if(inmap["-o"].size() > 2) throw string("\n More than one entry for " + inmap["-o"][0] + " option \n");
			fileout=inmap["-o"][1];
			fout=ofstream(fileout.c_str(),ios::out);
		}
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}

	gFin=&fin;
	gFout=&fout;
}
TrjRead::~TrjRead() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */

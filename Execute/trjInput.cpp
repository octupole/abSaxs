/*
 * trjInput.cpp
 *
 *  Created on: May 22, 2012
 *      Author: marchi
 */

#include "trjInput.h"

namespace trj {
trjInput::trjInput(int ntot,char ** v) {
	vector<string> in;

	inmap["-o"]=in;
	inmap["-i"]=in;
	inmap["-supcell"]=in;
	inmap["-grid"]=in;
	inmap["-dens"]=in;
	inmap["-qhist"]=in;
	inmap["-pdb"]=in;

	map<string,vector<string> >::iterator it=inmap.begin();
	for(int n=0;it!=inmap.end();++it,n++){
		Usage[n]=" ";
	}
	int M{0};
	Usage[M++]="\t -o fileout \n";
	Usage[M++]="\t -i filein \n";
	Usage[M++]="\t -supcell <double l>\n";
	Usage[M++]="\t -grid <double nx, [double ny, double nz] \n";
	Usage[M++]="\t -qhist <float dq=0.05> <float qcut=4>\n"
			"\t\tHistogram parameters as in -saxs, but for the difference calculation.\n";

	Usage[M++]="\t -dens <string select='R'> <int order=4> <int avg=2>\n"
			"\t\t Compute electron density instead of SAXS. select=[i++] R compute on R-space\n "
			"\t\t Q compute it from Q-space; order is the Lagrangian order; avg is how many bins\n"
			"\t\t are averaged\n";
	Usage[M++]="\t -pdb <string filename>\n"
			"\t\tInput a PDB file containing the entire system investigated. This input file \n "
			"\t\tis required when computing the SAXS intensity of the solvated protein or of the buffer. From this\n"
			"\t\tfile the topology of the system is extracted.\n";

	int n=1;
	string key;
	for(;n<ntot;){
		string tmp0(v[n]);
		if(v[n][0] =='-'){
			key.assign(v[n]);
			if(inmap.find(key) != inmap.end()){
				if(inmap[key].empty()) inmap[key].push_back(tmp0);
			} else{
				unknownC.push_back(key);
				inmap[key].push_back(tmp0);
			}
		}
		else{
			inmap[key].push_back(tmp0);
		}
		n++;
	}

}

vector<string> trjInput::getUsage(){
	vector<string> use;
	map<int,string>::iterator it=Usage.begin();
	for(;it!=Usage.end();++it){
		use.push_back(it->second);
	}
	return use;
}
trjInput::~trjInput() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */

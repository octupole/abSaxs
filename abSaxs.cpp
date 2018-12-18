/*
 * abSaxs.cpp
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */
#include <iostream>
#include "TrjRead.h"
#include "ExecabSaxs.h"
#include <mpi.h>
#include "fftw3-mpi.h"
#include "Topol.h"
#include "Atoms.h"

#include "HeaderTrj.h"
// #include "ResidueCM.h"
#include "TopolPDB.h"


using std::cout;
using std::endl;
using AtomsD=Atoms<double>;
using MAtoms=AtomsD;

using namespace Topol_NS;
MAtoms * atm{nullptr};

TopolPDB topPDB;

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
    fftw_mpi_init();

	/**
	 * read input
	 */

	ExecabSaxs * MyRun{nullptr};
	trj::TrjRead MyIn(argc,argv);
	MyIn.Input();

	Topol MyTop;

	vector<string> data;
	if(MyIn.gFpdb()){
		// read pdb file to construct topology
		for(string str;getline(*MyIn.gFpdb(),str);){
			data.push_back(str);
		}
		topPDB(data);
		MyTop(topPDB,false);
		atm=new Atoms<double>(MyIn.WhichDiffusion);
		int natoms=MyTop.Size();
		atm->setDim(natoms);
		atm->setTopol(MyTop);
		atm->initLists(topPDB);
	}


	cout << "Here "<<endl;
	MyRun=new ExecabSaxs(MyIn,MyTop);
	MyRun->Run_abSaxs(atm);
    MPI_Finalize();
    return 0;
}



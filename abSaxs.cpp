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

using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
    fftw_mpi_init();

	/**
	 * read input
	 */

	ExecabSaxs * MyRun{nullptr};
	trj::TrjRead MyIn(argc,argv);
	MyIn.Input();
	MyRun=new ExecabSaxs(MyIn);
	MyRun->Run_abSaxs();
    MPI_Finalize();
    return 0;
}



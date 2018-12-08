/*
 * abSaxs.cpp
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */
#include <iostream>
#include "TrjRead.h"
#include "ExecabSaxs.h"
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	/**
	 * read input
	 */
	ExecabSaxs * MyRun{nullptr};
	trj::TrjRead MyIn(argc,argv);
	MyIn.Input();
	MyRun=new ExecabSaxs(MyIn);
	MyRun->Run_abSaxs();

}



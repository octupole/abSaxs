/*
 * abSaxs.cpp
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */
#include <iostream>
#include "TrjRead.h"
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	/**
	 * read input
	 */
	trj::TrjRead MyIn(argc,argv);
	MyIn.Input();

  cout << "Ula " <<endl;
}



/*
 * TrjRead.h
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */

#ifndef SRC_TRJREAD_H_
#define SRC_TRJREAD_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>

#include "trjInput.h"
#include "myEnums.hpp"
#include "MyUtilClass.h"
#include "DensMode.h"
#include "CenterMass.h"



using namespace std;

namespace trj {
template <typename T>
class Streams{
	T * myStream{nullptr};
public:
	Streams(){};
	Streams(T * y):myStream{&y}{};
	 Streams & operator=(T * y ){
		myStream=y;
		return *this;
	}
	T * operator()(){
		return myStream;
	}

};

template <typename T>
class Values{
	T * const myValue;
public:
	Values(T & y): myValue{&y}{};
	T & operator()(){return *myValue;}
};

/** \brief The class extends trjInput by providing methods to read the on line options
 *         and provide the stored data to calling classes.
 *
 */
using Dvect=DVECT::DDvect<double>;
class TrjRead: public trjInput {
	int MyOrder{1};
	int MyDensAvg{0};
	size_t BoxMultiply{1};
	double Mycut{4.0};
	double Myd{0.05};
	string filein,fileout{"abSaxs.out"},filepdb,filedens{"abSaxs.vts"};
	double SuperCellSide{1.0};
	unsigned int nnx{128},nny{128},nnz{128};
	ifstream * finx{nullptr};
	ifstream * fpdb{nullptr};
	ofstream * foutx{nullptr};
	DensMode ModeCompute;

public:
	Enums::Compute WhatToDo{Enums::SAXS};

	Values<string> gfilein{filein};
	Values<string> gfileout{fileout};
	Values<string> gfiledens{filedens};
	Values<double> gSupCell{SuperCellSide};
	Streams<ifstream> gFin;
	Streams<ofstream> gFout;
	Streams<ifstream> gFpdb;
	Values<double> gMyCut{Mycut};
	Values<double> gMyd{Myd};
	CenterMass_t WhichDiffusion{diffk};

	Values<unsigned int> gnnx{nnx};
	Values<unsigned int> gnny{nny};
	Values<unsigned int> gnnz{nny};
	Values<int> gMyOrder{MyOrder};
	Values<int> gMyDensAvg{MyDensAvg};
	Values<DensMode> gModeCompute{ModeCompute};

	TrjRead(int nv,char ** v);
	void Input();
	virtual ~TrjRead();
};

} /* namespace trj */

#endif /* SRC_TRJREAD_H_ */

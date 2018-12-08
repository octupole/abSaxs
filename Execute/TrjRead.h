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
	string filein,fileout{"abSaxs.out"};
	ifstream * finx{nullptr};
	ofstream * foutx{nullptr};
public:
	Values<string> gfilein{filein};
	Values<string> gfileout{fileout};
	Streams<ifstream> gFin;
	Streams<ofstream> gFout;

	TrjRead(int nv,char ** v);
	void Input();
	virtual ~TrjRead();
};

} /* namespace trj */

#endif /* SRC_TRJREAD_H_ */

/*
 * ExecabSaxs.h
 *
 *  Created on: Dec 8, 2018
 *      Author: marchi
 */

#ifndef EXECUTE_EXECABSAXS_H_
#define EXECUTE_EXECABSAXS_H_

#include "TrjRead.h"
#include "MyUtilClass.h"
#include "SGSmooth.h"
#include <cstdlib>
#include <vector>
#include <sstream>

using std::vector;

class ExecabSaxs {
	ifstream * ExIn;
	ofstream * ExOut;
public:
	ExecabSaxs(trj::TrjRead &);
	ExecabSaxs()=delete;
	void Run_abSaxs();
	virtual ~ExecabSaxs();
};

#endif /* EXECUTE_EXECABSAXS_H_ */

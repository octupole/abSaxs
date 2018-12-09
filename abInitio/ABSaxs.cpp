/*
 * ABSaxs.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#include "ABSaxs.h"
namespace abinit{
ABSaxs::ABSaxs() {
}

void ABSaxs::setUp(SaxsData * exp){
	double Rg=exp->getRg();
	double Rd=exp->Rd();
	metric(Rd);
	cout << CO[XX][XX] <<endl;
}
void ABSaxs::Run(){

}
ABSaxs::~ABSaxs() {
	// TODO Auto-generated destructor stub
}
}

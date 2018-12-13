/*
 * Funktionell.cpp
 *
 *  Created on: Dec 13, 2018
 *      Author: marchi
 */

#include "Funktionell.h"

namespace Funkll {

Funktionell::Funktionell(RhoSaxs & ro_out, RhoSaxs & ro_in){
	CO=ro_out.getCO();
	OC=ro_out.getOC();
	co=ro_in.getCO();
	oc=ro_in.getOC();
	nx=ro_out.getnnx();
	ny=ro_out.getnny();
	nz=ro_out.getnnz();
	nzp=nz/2+1;
	Nx=ro_in.getnnx();
	Ny=ro_in.getnny();
	Nz=ro_in.getnnz();
}
double Funktionell::Deviate(array3<Complex> & InC){
	double energy{0};

	return energy;
}
Funktionell::~Funktionell() {
}

} /* namespace Funkll */

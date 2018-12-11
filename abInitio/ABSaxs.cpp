/*
 * ABSaxs.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#include "ABSaxs.h"
namespace abinit{

ABSaxs::ABSaxs(uint nx, uint ny, uint nz, double SupCell=3.0): grid_a{nx,ny,nz},
		SuperCell{SupCell} {
		}

void ABSaxs::setUp(SaxsData * exp){
	double Rg=exp->getRg();
	double Rd=exp->Rd();
	cout << "Found Rg = "<< Rg << " " <<"Rd = "<<Rd<<endl;
	metric(Rd);
	std::function<uint(uint)> cell=[this](uint n){double tmp=(double) n/SuperCell;
		uint tmp0=(uint) tmp;return tmp0;};
	grid_b[XX]=cell(grid_a[XX]);
	grid_b[YY]=cell(grid_a[YY]);
	grid_b[ZZ]=cell(grid_a[ZZ]);

	CO[XX][XX]=(co[XX][XX]/(double) grid_b[XX])*(double) grid_a[XX];
	CO[YY][YY]=(co[YY][YY]/(double) grid_b[YY])*(double) grid_a[YY];
	CO[ZZ][ZZ]=(co[ZZ][ZZ]/(double) grid_b[ZZ])*(double) grid_a[ZZ];

}
void ABSaxs::Run(){
	Rho_in.Allocate(grid_a[XX],grid_a[YY],grid_a[ZZ]);
	Rho_s.Allocate(grid_b[XX],grid_b[YY],grid_b[ZZ]);
}
ABSaxs::~ABSaxs() {
	// TODO Auto-generated destructor stub
}
}

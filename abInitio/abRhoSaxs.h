/*
 * RhoSaxs.h
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#ifndef SRC_ABINIT_RHOSAXS_H_
#define SRC_ABINIT_RHOSAXS_H_
//#ifdef HAVE_VTK
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkFloatArray.h"
#include "vtkHedgeHog.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkStructuredGrid.h"
#include "vtkXMLStructuredGridWriter.h"
//#endif

#include <map>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>

#include "myEnums.hpp"
#include "Grid.h"
#include "MyUtilClass.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "Finalize.h"
using namespace MATRIX;
using namespace DVECT;

using std::map;

using namespace Enums;
/** \brief Derived class of Grid<1>
 *         A tridimensional grid storing the atomic density.
 */
namespace abInitioRho{
const int NMult=1;
class RhoSaxs: public Grid<1> {
protected:
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	static bool firstTime;
	Matrix CO,OC;
	double dx{0},dy{0},dz{0};

	size_t Nx{0},Ny{0},Nz{0};
	bool __check(){
		bool notOk{false};
		for(int n{0};n<DIM;n++)
			if(CO[n]) notOk=true;
		return notOk;
	}
	void __setDx();
public:
	RhoSaxs(Matrix & co):CO{co}{
		OC=CO.Inversion();__setDx();};
	RhoSaxs()=delete;
	RhoSaxs(const RhoSaxs & y): Grid<1>::Grid<1>(y){
		CO=y.CO;
		OC=y.OC;
		Nx=y.Nx;
		Ny=y.Ny;
		Nz=y.Nz;
		__setDx();
	}
	RhoSaxs(size_t nx,size_t ny,size_t nz,Matrix & co):CO{co},Grid<1>::Grid<1>(nx,ny,nz)
			{OC=CO.Inversion();__setDx();};
	void setCO(Matrix CO){this->CO=CO;this->OC=this->CO.Inversion();__setDx();}
	void PartialCopy(RhoSaxs &);
	void MakeAvg();
	void setNx(size_t Mx,size_t My, size_t Mz){
		Nx=Mx;Ny=My;Nz=Mz;
	}
	Matrix & getCO(){return CO;}
	Matrix & getOC(){return OC;}
	virtual RhoSaxs & operator=(const double y){
		this->Grid<1>::operator =(y);
		return *this;
	};
	virtual RhoSaxs & operator=(const RhoSaxs & y){
		this->Grid<1>::operator =(y);
		return *this;
	};
	void initDensity(double);
	void initBarrier(double,double,array3<double> &);
	void Density();
	void WriteIt();
	void copyIn(array3<double> &);
	void copyOut(array3<double> &);
	virtual ~RhoSaxs(){};
	friend std::ofstream & operator<<(std::ofstream &,RhoSaxs & );
};
}
#endif /* SRC_ABINIT_RHOSAXS_H_ */

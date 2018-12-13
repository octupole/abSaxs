/*
 * RhoSaxs.h
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#ifndef SRC_RHOSAXS_H_
#define SRC_RHOSAXS_H_
#ifdef HAVE_VTK
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
#endif

#include <map>
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
const int NMult=1;

using namespace Enums;
/** \brief Derived class of Grid<1>
 *         A tridimensional grid storing the atomic density.
 */
namespace abinit{
class RhoSaxs: public Grid<1> {
protected:
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	static bool firstTime;
public:
	RhoSaxs(){};
	RhoSaxs(const RhoSaxs & y): Grid<1>::Grid<1>(y){}
	RhoSaxs(size_t nx,size_t ny,size_t nz):Grid<1>::Grid<1>(nx,ny,nz) {};
	void MakeAvg();
	virtual RhoSaxs & operator=(const double y){
		this->Grid<1>::operator =(y);
		return *this;
	};
	virtual RhoSaxs & operator=(const RhoSaxs & y){
		this->Grid<1>::operator =(y);
		return *this;
	};

	void Density();
	virtual void Write();
	virtual ~RhoSaxs(){};
};
}
#endif /* SRC_RHOSAXS_H_ */

/*
 * Prcfft3d.h
 *
 *  Created on: Dec 7, 2018
 *      Author: marchi
 */

#ifndef FFTW3PP_PRCFFT3D_H_
#define FFTW3PP_PRCFFT3D_H_

#include "fftw3-mpi.h"

#include "Array.h"
#include "Pfftw.h"

namespace Pfftwpp{
class Prcfft3d: public Pfftw {
	  unsigned int nx;
	  unsigned int ny;
	  unsigned int nz;
public:
	Prcfft3d(unsigned int nx, unsigned int ny, unsigned int nz, const Array::array3<double>& in,
			const Array::array3<Complex>& out)
    : Pfftw(in.Size(),nx*ny*nz), nx(nx), ny(ny), nz(nz) {
		Setup(in,out);
	}

	fftw_plan Plan(Complex *in, Complex *out) {

		return  fftw_mpi_plan_dft_r2c_3d(nx,ny,nz,(double *) in,(fftw_complex *) out,MPI_COMM_WORLD,
				effort);
	}

	void Execute(Complex *in, Complex *out, bool shift=false) {
		fftw_execute(plan);
	}

};
}
#endif /* FFTW3PP_PRCFFT3D_H_ */

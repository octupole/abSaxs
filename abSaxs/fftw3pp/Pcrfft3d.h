/*
 * Pcrfft3d.h
 *
 *  Created on: Dec 7, 2018
 *      Author: marchi
 */

#ifndef FFTW3PP_PCRFFT3D_H_
#define FFTW3PP_PCRFFT3D_H_

#include "fftw3-mpi.h"

#include "Pfftw.h"
#include "Array.h"
namespace fftwpp {

class Pcrfft3d: public Pfftw {
	  unsigned int nx;
	  unsigned int ny;
	  unsigned int nz;
public:
	Pcrfft3d(unsigned int nx, unsigned int ny, unsigned int nz, const Array::array3<Complex>& in,
			const Array::array3<double>& out):
		Pfftw(out.Size(),nx*ny*nz), nx(nx), ny(ny), nz(nz) {Setup(in,out);}

	fftw_plan Plan(Complex *in, Complex *out) {
		fftw_plan tmp=fftw_mpi_plan_dft_c2r_3d(nx,ny,nz,(fftw_complex *) in,(double *) out,
		                                MPI_COMM_WORLD,effort);
	return tmp;
}

void Execute(Complex *in, Complex *out, bool shift=false) {
	fftw_execute(plan);
}

};

} /* namespace fftwpp */

#endif /* FFTW3PP_PCRFFT3D_H_ */

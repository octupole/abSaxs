/*
 * Pfftw.h
 *
 *  Created on: Dec 7, 2018
 *      Author: marchi
 */

#ifndef FFTW3PP_PFFTW_H_
#define FFTW3PP_PFFTW_H_
#include <fstream>
#include <iostream>

#include "fftw3-mpi.h"
#include <complex>

#include "Array.h"
using Complex=std::complex<double>;

inline void fftwpp_export_wisdom(void (*emitter)(char c, std::ofstream& s),
                                 std::ofstream& s)
{
  fftw_export_wisdom((void (*) (char, void *)) emitter,(void *) &s);
}

inline int fftwpp_import_wisdom(int (*g)(std::ifstream& s), std::ifstream &s)
{
  return fftw_import_wisdom((int (*) (void *)) g,(void *) &s);
}

inline void PutWisdom(char c, std::ofstream& s) {s.put(c);}
inline int GetWisdom(std::ifstream& s) {return s.get();}



namespace fftwpp {

class Pfftw {
protected:
	unsigned int doubles; // number of double words in dataset

	double norm;
	int rank{0},size{0};
	fftw_plan plan{nullptr};

	static std::ifstream ifWisdom;
	static std::ofstream ofWisdom;
	static bool Wise;
	static const double twopi;

public:
	static unsigned int effort;
	static const char *WisdomName;
	Pfftw(unsigned int doubles,unsigned int n):
		doubles{doubles},norm(1.0/(double) n){
			if(!size){
				MPI_Comm_size(MPI_COMM_WORLD,&size);
				MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			}
		}

	virtual ~Pfftw() {if(plan) fftw_destroy_plan(plan);}

	virtual fftw_plan Plan(Complex *in, Complex *out)=0;

	virtual void fftNormalized(Complex *in, Complex *out=NULL) {
		Execute(in,out);
		Normalize(out);
	}
	virtual void Execute(Complex *in, Complex *out, bool=false) {
		fftw_execute_dft(plan,(fftw_complex *) in,(fftw_complex *) out);
	}

	void Setup(Complex *in, Complex *out=NULL) {
		if(!Wise)LoadWisdom();

		plan=Plan(in,out);
		if(!plan) {
			std::cerr << "Unable to construct FFTW plan" << std::endl;
			exit(1);
		}
		SaveWisdom();

	}

	void Setup(Complex *in, double *out) {Setup(in,(Complex *) out);}

	void Setup(double *in, Complex *out=NULL) {Setup((Complex *) in,out);}

	void LoadWisdom() {
		if(!rank){
			ifWisdom.open(WisdomName);
			fftwpp_import_wisdom(GetWisdom,ifWisdom);
			ifWisdom.close();
		}
		fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);
		Wise=true;
	}

	void SaveWisdom() {
		fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
		if(rank == 0){
			ofWisdom.open(WisdomName);
			fftwpp_export_wisdom(PutWisdom,ofWisdom);
			ofWisdom.close();
		}
	}

	void Normalize(Complex *out) {
		unsigned int stop=(doubles+1)/2;
		for(unsigned int i=0; i < stop; i++) out[i] *= norm;
	}

	void Normalize(double *out) {
		for(unsigned int i=0; i < doubles; i++) out[i] *= norm;
	}


	void fft(Complex *in, Complex *out=NULL) {
		Execute(in,out);
	}

	void fft(double *in, Complex *out=NULL) {
		fft((Complex *) in,out);
	}

	void fft(Complex *in, double *out) {
		fft(in,(Complex *) out);
	}
	void fftNormalized(Complex *in, double *out) {
		Execute(in,(Complex *) out);
		Normalize(out);
	}

	void fftNormalized(double *in, Complex *out) {
		fftNormalized((Complex *) in,out);
	}

};

} /* namespace fftwpp */

#endif /* FFTW3PP_PFFTW_H_ */

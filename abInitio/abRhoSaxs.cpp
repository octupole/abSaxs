/*
 * RhoSaxs.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#include "abRhoSaxs.h"
namespace abInitioRho{
bool RhoSaxs::firstTime=true;
void RhoSaxs::Density(){
}
void RhoSaxs::initDensity(double Rd /* Radius of the sphere not Rg!*/){
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution (0.2,1.0);

	try{
		if(!__check()) throw string("\nCO matrix not initialized.\n");
		if(this->Size() == 0) throw string("\n Electron density: array not allocated.\n");
	}catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	Dvect P=CO*Dvect{0.5,0.5,0.5};
	/* Assume Cell axis are orthogonal */
	double dx=CO[XX][XX]/(double) this->getnnx();
	double dy=CO[YY][YY]/(double) this->getnny();
	double dz=CO[ZZ][ZZ]/(double) this->getnnz();
	double dens{0.6},Rd2=Rd*Rd;
	for(int o{0};o< this->getnnx();o++){
		double xc=dx*(double)o;
		for(int p{0};p<this->getnny();p++){
			double yc=dy*(double)p;
			for(int q{0};q< this->getnnz();q++){
				double zc=dz*(double)q;
				Dvect grid{xc,yc,zc};

				if((grid-P).Norm2() <=Rd2)
				(*this)[0][o][p][q]=distribution(generator);

			}
		}
	}
}
void RhoSaxs::initBarrier( double Pot, double Rd,array3<double> & x){
	try{
		if(!__check()) throw string("\nCO matrix not initialized.\n");
		if(this->Size() == 0) throw string("\n Electron density: array not allocated.\n");
	}catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	Dvect P=CO*Dvect{0.5,0.5,0.5};
	/* Assume Cell axis are orthogonal */
	double dx=CO[XX][XX]/(double) this->getnnx();
	double dy=CO[YY][YY]/(double) this->getnny();
	double dz=CO[ZZ][ZZ]/(double) this->getnnz();
	double dens{0.2},Rd2=(Rd+6)*(Rd+6);
	for(int o{0};o< this->getnnx();o++){
		double xc=dx*(double)o;
		for(int p{0};p<this->getnny();p++){
			double yc=dy*(double)p;
			for(int q{0};q< this->getnnz();q++){
				double zc=dz*(double)q;
				Dvect grid{xc,yc,zc};

				if((grid-P).Norm2() <=Rd2)
					x[o][p][q]=Pot;

			}
		}
	}

}

void RhoSaxs::PartialCopy(RhoSaxs & rho0){
	if(this->Size()>rho0.Size()){
		for(int o{0};o<rho0.getnnx();o++){
			for(int p{0};p<rho0.getnny();p++){
				for(int q{0};q< rho0.getnnz();q++){
					(*this)[0][o][p][q]=rho0[0][o][p][q];
				}
			}
		}
	} else{
		for(int o{0};o<this->getnnx();o++){
			for(int p{0};p<this->getnny();p++){
				for(int q{0};q< this->getnnz();q++){
					(*this)[0][o][p][q]=rho0[0][o][p][q];
				}
			}
		}
	}
}
void RhoSaxs::__setDx(){
	dx=CO[XX][XX]/(double) this->getnnx();
	dy=CO[YY][YY]/(double) this->getnny();
	dz=CO[ZZ][ZZ]/(double) this->getnnz();
}

void RhoSaxs::copyIn(array3<double> & ro){
	unsigned int nx=this->getnnx();
	unsigned int ny=this->getnny();
	unsigned int nz=this->getnnz();
	for(size_t o{0};o<nx;o++)
		for(size_t p{0};p<ny;p++)
			for(size_t q{0};q<nz;q++){
				(*this)[0][o][p][q]=ro[o][p][q];
			}
}
void RhoSaxs::copyOut(array3<double> & ro){
	unsigned int nx=this->getnnx();
	unsigned int ny=this->getnny();
	unsigned int nz=this->getnnz();
	for(size_t o{0};o<nx;o++)
		for(size_t p{0};p<ny;p++)
			for(size_t q{0};q<nz;q++){
				ro[o][p][q]=(*this)[0][o][p][q];
			}
}
void RhoSaxs::WriteIt(){
	double dvol=static_cast<float>(dx*dy*dz);
	unsigned int mx=this->getnnx();
	unsigned int my=this->getnny();
	unsigned int mz=this->getnnz();
	Nx=Nx==0?mx:Nx;
	Ny=Ny==0?my:Ny;
	Nz=Nz==0?mz:Nz;
	float x[3], a, v[3], rMin = 0.5, rMax = 1., deltaRad, deltaZ;
	int dims[3]={static_cast<int>(Nx),static_cast<int>(Ny),static_cast<int>(Nz)};
	  // Create the structured grid.
	vtkStructuredGrid *sgrid = vtkStructuredGrid::New();
 	sgrid->SetDimensions(dims);

	  // Create points.
	vtkPoints *points = vtkPoints::New();
	points->Allocate(dims[0]*dims[1]*dims[2]);

	  // Create scalars. No need to allocate size, as the array will expand automatically as we add new
	  // values with InsertNextValue.
	vtkFloatArray *scalars = vtkFloatArray::New();
	scalars->Allocate(dims[0]*dims[1]*dims[2]);
	scalars->SetName("density");
	for (size_t k{0}; k < Nz; k++) {
		x[ZZ]=dz*k;
		int kOffset = k*Nx*Ny;
	    for (size_t j{0}; j<Ny; j++) {
	    	int jOffset = j * Nx;
			x[YY]=dy*j;
	    	for (size_t i{0}; i<Nx; i++) {
	    		int offset = i + jOffset + kOffset;
	    		x[XX]=dx*i;

	    		float Val=(*this)[0][i][j][k]/dvol;
	    		points->InsertPoint(offset,x);
	    		scalars->InsertValue(offset,Val);
	    	}
	    }
	}
	sgrid->SetPoints(points);
	points->Delete();
	sgrid->GetPointData()->SetScalars(scalars);
	scalars->Delete();
	vtkXMLStructuredGridWriter *Writer = vtkXMLStructuredGridWriter::New();
	Writer->SetInputData(sgrid);
	Writer->SetFileName("RhoSaxs.vts");
	Writer->Write();

	// Delete objects.
	Writer->Delete();
	sgrid->Delete();

//#else
//	for(size_t i{0};i<nx;i++)
//		for(size_t j{0};j<ny;j++)
//			for(size_t k{0};k<nz;k++){
//				fout << dx*i<< ", " << dy*j<< ", " << dz*k<< ", " << I_r0[i][j][k] << endl;
//			}
//#endif

}

std::ofstream & operator<<(std::ofstream & fout, RhoSaxs & y){

	return fout;
}
}

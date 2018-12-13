/*
 * ABSaxs.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#include "ABSaxs.h"
namespace abinit{

ABSaxs::ABSaxs(uint nx, uint ny, uint nz, double SupCell): grid_b{nx,ny,nz},
		nx{nx},ny{ny},nz{nz},nzp{nz/2+1},SuperCell{SupCell} {
		}

void ABSaxs::setUp(SaxsData * exp){
	double Rg=exp->getRg();
	double Rd=exp->Rd();
	metric(Rd);
	cout << " Rg = " << Rg << " Rd ="<<Rd<<endl;
	std::function<uint(uint)> cell=[this](uint n){double tmp=(double) n/SuperCell;
		uint tmp0=(uint) tmp;return tmp0;};
	grid_a[XX]=cell(grid_b[XX]);
	grid_a[YY]=cell(grid_b[YY]);
	grid_a[ZZ]=cell(grid_b[ZZ]);
	Nx=grid_a[XX];
	Ny=grid_a[YY];
	Nz=grid_a[ZZ];


	CO[XX][XX]=(co[XX][XX]/(double) grid_a[XX])*(double) grid_b[XX];
	CO[YY][YY]=(co[YY][YY]/(double) grid_a[YY])*(double) grid_b[YY];
	CO[ZZ][ZZ]=(co[ZZ][ZZ]/(double) grid_a[ZZ])*(double) grid_b[ZZ];
	OC=CO.Inversion();
	oc=co.Inversion();
	dx=CO[XX][XX]/static_cast<double>(nx);
	dy=CO[YY][YY]/static_cast<double>(ny);
	dz=CO[ZZ][ZZ]/static_cast<double>(nz);
	dvol=static_cast<float>(dx*dy*dz);

	Rho_in=new RhoSaxs(grid_a[XX],grid_a[YY],grid_a[ZZ],co);
	Rho_s=new RhoSaxs(grid_b[XX],grid_b[YY],grid_b[ZZ],CO);

	Rho_in->initDensity(Rd);
	Rho_s->CopySmallerRho(*Rho_in);

	ptrdiff_t alloc_local, local_n0, local_0_start;

     alloc_local = fftw_mpi_local_size_3d(nx, ny, nzp, MPI_COMM_WORLD,
                                            &local_n0, &local_0_start);

 	F_k.Allocate(grid_b[XX],grid_b[YY],nzp);
	F_r.Allocate(grid_b[XX],grid_b[YY],nzp*2);
	I_k.Allocate(grid_b[XX],grid_b[YY],nzp);
	cout << "Run with Rho_inner ="<< Nx<<" "<< Ny<<" "<< Nz<<" Resolution "<< co[XX][XX]/(float) Nx<<endl;
	cout << "Run with Rho_outer ="<< nx<<" "<< ny<<" "<< nz<<" "<<endl;
}
void ABSaxs::Run(){
	Minimize();

}
array3<Complex> ABSaxs::Modulus(array3<Complex> & ro_k){
	array3<Complex> ro_k1(nx,ny,nzp);
#pragma omp parallel for collapse(3)
	for(auto i=0;i<nx;i++){
		for(auto j=0;j<ny;j++){
			for(auto k=0;k<nzp;k++){
				ro_k1[i][j][k]=ro_k[i][j][k]*conj(ro_k[i][j][k]);
			}
		}
	}
	return ro_k1;
}

void ABSaxs::Minimize(){
	/*
	 * Arrays are destroyed while performing plan!!!
	 */
	Pfftwpp::Pcrfft3d Backward3(nx,ny,nz,F_k,F_r);
	Pfftwpp::Prcfft3d  Forward3(nx,ny,nz,F_r,F_k);

	Rho_s->copyOut(F_r);
	F_r*=dvol;

	Forward3.fft(F_r,F_k);
	I_k=Modulus(F_k);
	__qhistogram();
}

void ABSaxs::__qhistogram(){
	Matrix oc{OC};
	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};

	double mw1,mw2,mw3,mw;
	Dvect fx{(double)nfx-1,(double)nfy-1,(double)nfz-1};
	double m_qcut{this->qcut},m_dq{this->dq};
	vector<double> mydq0={2.0*M_PI*oc[XX][XX],2.0*M_PI*oc[YY][YY],2.0*M_PI*oc[ZZ][ZZ],this->dq};
	vector<double> mycut0={2.0*M_PI*oc[XX][XX]*fx[XX],2.0*M_PI*oc[YY][YY]*fx[YY],2.0*M_PI*oc[ZZ][ZZ]*fx[ZZ],this->qcut};

	double dq=1.1*(*std::max_element(mydq0.begin(),mydq0.end()));
	double qcut=*std::min_element(mycut0.begin(),mycut0.end());

	try{
		int ntry{0};
		bool notok=false;
		stringstream ss;
		string msg("");
		if(dq != m_dq){
			ss<< "from " <<m_dq << " to "<<dq ;
			msg+="    Histogram was constructed with a larger bin " + ss.str();
			notok=true;
			ntry++;
		}
		if(qcut != m_qcut) {
			ss.str(string());
			ss<< "from " <<m_qcut << " to "<<qcut ;
			if(ntry) msg+="\n";
			msg+="    Cutoff was too large. Changed to maximum allowed " + ss.str();
			notok=true;
			this->qcut=qcut;
			ntry++;
		}
		if(notok) throw msg;
	}catch(const string & s){
		cout << "\n***"+string(80,'-')+"***\n";
		cout << s <<endl;
		cout << "***"+string(80,'-')+"***\n";
	}
	map<unsigned int,double> qdf;
	map<uint,uint> qdn;
	cout << dq << " " <<endl;

	for(auto i=0;i<nx;i++){
		int ia=(i<nfx)?i : i-nx;
		size_t ib=i==0?0:nx-i;
		for(auto j=0;j<ny;j++){
			int ja=(j<nfy)?j : j-ny;
			size_t jb=j==0?0:ny-j;
			for(auto k=0;k<nzp;k++){
				int ka=(k<nfz)?k : k-nz;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw1=2.0*M_PI*mw1;
				mw2=2.0*M_PI*mw2;
				mw3=2.0*M_PI*mw3;
				mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
				if(mw<qcut){
					Complex v0=I_k[i][j][k];
					if(k != 0 && k != nzp-1){
						Complex vt1=I_k[ib][jb][k];
						v0=0.5*(v0+vt1);
					}
					int h0=static_cast<int>(mw/dq);
					int h1=h0+1;

					qdf[h0]+=v0.real();
					qdn[h0]++;
					if(h0 != 0){
						qdf[h1]+=v0.real();
						qdn[h1]++;
					}
				}
			}
		}
	}
	for(auto it=qdf.begin();it!=qdf.end();it++){
		auto h0=it->first;
		qdf[h0]/=qdn[h0];
		cout << h0*dq << " " << qdf[h0]<<endl;
	}


}
ABSaxs::~ABSaxs() {
	// TODO Auto-generated destructor stub
}
}

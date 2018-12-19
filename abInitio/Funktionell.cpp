/*
 * Funktionell.cpp
 *
 *  Created on: Dec 13, 2018
 *      Author: marchi
 */

#include "Funktionell.h"
const double INEN{100.0};
namespace Funkll {
struct pickPar{
	pickPar(double &Par, double & scl,map<size_t,double> & Ie,map<size_t,double> & Ic){
		double energy{0};
		auto i0=Ie.begin();
		double sc=Ic[i0->first]/Ie[i0->first];
		for(auto it{Ie.begin()};it!=Ie.end();it++){
			auto h0=it->first;
			Ie[h0]*=sc;
			double tmp=Ie[h0]-Ic[h0];
			double wei=1.0;
			energy+=tmp*tmp;
		}
		scl=1.0;
		Par=INEN/energy;
	}

};


Funktionell::Funktionell(abInitioRho::RhoSaxs * ro_out, abInitioRho::RhoSaxs * ro_in,SaxsData * exp){
	CO=ro_out->getCO();
	OC=ro_out->getOC();
	co=ro_in->getCO();
	oc=ro_in->getOC();
	nx=ro_out->getnnx();
	ny=ro_out->getnny();
	nz=ro_out->getnnz();
	nzp=nz/2+1;
	Nx=ro_in->getnnx();
	Ny=ro_in->getnny();
	Nz=ro_in->getnnz();
	Rd=co[XX][XX]*2.0/3.5;
	this->setUpFirst(exp);
}
void Funktionell::setUpFirst(SaxsData * exp){
	double m_qcut{this->qcut},m_dq{this->dq};
	vector<double> mydq0={2.0*M_PI*OC[XX][XX],2.0*M_PI*OC[YY][YY],2.0*M_PI*OC[ZZ][ZZ],this->dq};
	double dq=1.1*(*std::max_element(mydq0.begin(),mydq0.end()));
	this->dq=dq;

	vector<std::pair<double,double>> & tmp=exp->gIq_exp();
	map<size_t,double> qft;
	map<size_t,size_t> nqft;
	for(size_t o{0};o<tmp.size();o++){
		int h0=static_cast<int>(tmp[o].first/dq);
		qft[h0]+=tmp[o].second;
		nqft[h0]++;
	}
	qmin=qft.begin()->first*dq;
	for(auto it{qft.begin()};it!=qft.end();it++){
		const int h0=it->first;
		Iq_exp[h0]=qft[h0]/static_cast<double>(nqft[h0]);
		qcut=h0*dq+0.5*dq;
	}

	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};
	double mw1,mw2,mw3,mw;

	for(auto i=0;i<nx;i++){
		int ia=(i<nfx)?i : i-nx;
		size_t ib=i==0?0:nx-i;
		for(auto j=0;j<ny;j++){
			int ja=(j<nfy)?j : j-ny;
			size_t jb=j==0?0:ny-j;
			for(auto k=0;k<nzp;k++){
				int ka=(k<nfz)?k : k-nz;
				mw1=OC[XX][XX]*ia+OC[XX][YY]*ja+OC[XX][ZZ]*ka;
				mw2=OC[YY][XX]*ia+OC[YY][YY]*ja+OC[YY][ZZ]*ka;
				mw3=OC[ZZ][XX]*ia+OC[ZZ][YY]*ja+OC[ZZ][ZZ]*ka;
				mw1=2.0*M_PI*mw1;
				mw2=2.0*M_PI*mw2;
				mw3=2.0*M_PI*mw3;
				mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
				int h0=static_cast<int>(mw/dq);
				if(mw<qcut && mw >=qmin){
					mapIdx[h0].push_back(vector<int>{i,j,k});
				}
			}
		}
	}
}
double Funktionell::EnergyR(array3<double> & F_r,array3<double> & grad){
	double energy{0};
	Dvect P=co*Dvect{0.5,0.5,0.5};
	/* Assume Cell axis are orthogonal */
	double dx=co[XX][XX]/(double) Nx;
	double dy=co[YY][YY]/(double) Ny;
	double dz=co[ZZ][ZZ]/(double) Nz;
	double R2=Rd*Rd;
	for(int o{0};o< Nx;o++){
		double xc=dx*(double)o;
		for(int p{0};p<Ny;p++){
			double yc=dy*(double)p;
			for(int q{0};q<Nz;q++){
				double zc=dz*(double)q;
				Dvect grid{xc,yc,zc};
				if((grid-P).Norm2()>R2){
					double tmp{F_r[o][p][q]};
					energy+=Pot*tmp*tmp;
					grad[o][p][q]+=2.0*Pot*tmp;
				}
			}
		}
	}

	return energy;
}
double Funktionell::EnergyQ(double & scalePlot,double & Gscale,array3<Complex> & F_k){
	array3<Complex> Grad(nx,ny,nzp,sizeof(Complex));
	Iq_c.clear();
	auto InC=Modulus(F_k);
	Complex vt0{0,0};
	double energy{0};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		const size_t h0=it->first;
		vector<vector<int>> & vIdx=it->second;
		for(size_t o{0};o<vIdx.size();o++){
			int i=vIdx[o][XX];
			int j=vIdx[o][YY];
			int k=vIdx[o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				vt0=InC[i][j][k]+InC[ib][jb][k];
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());

	}
	static pickPar Once(Par,scalePlot,Iq_exp,Iq_c);
	myScale=scalePlot;
	Grad=Complex{0,0};
	energy=0;
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		auto tmp=(Iq_exp[h0]-Iq_c[h0]*scalePlot);
		double wei=1.0;
		energy+=tmp*tmp*Par*wei;
		double grad=-2.0*Par*wei*tmp*scalePlot/(double) mapIdx[h0].size();
		Gscale+=-2.0*Par*wei*tmp*Iq_c[h0];
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				Grad[i][j][k]+=grad*F_k[i][j][k];
				Grad[ib][jb][k]+=grad*F_k[ib][jb][k];
			} else{
				Grad[i][j][k]+=2.0*grad*F_k[i][j][k];
			}
		}
	}
	F_k=Grad;
	return energy;
}
double Funktionell::EnergyQQ(double & scalePlot,double & Gscale,array3<Complex> & F_k){
	array3<Complex> Grad(nx,ny,nzp,sizeof(Complex));
	Iq_c.clear();
	auto InC=Modulus(F_k);
	Complex vt0{0,0};
	double energy{0};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		const size_t h0=it->first;
		vector<vector<int>> & vIdx=it->second;
		for(size_t o{0};o<vIdx.size();o++){
			int i=vIdx[o][XX];
			int j=vIdx[o][YY];
			int k=vIdx[o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				vt0=InC[i][j][k]+InC[ib][jb][k];
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());
	}
	static pickPar Once(Par,scalePlot,Iq_exp,Iq_c);
	myScale=scalePlot;
	Grad=Complex{0,0};
	energy=0;
	double tmp,grad;
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		double wei=1;
//		double wei=1.0/Iq_exp[h0];
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				tmp=Iq_exp[h0]-InC[i][j][k].real();
				energy+=tmp*tmp*Par*wei;
				grad=-2.0*Par*wei*tmp;
				Grad[i][j][k]+=grad*F_k[i][j][k];

				tmp=Iq_exp[h0]-InC[ib][jb][k].real();
				energy+=tmp*tmp*Par*wei;
				grad=-2.0*Par*wei*tmp*scalePlot;
				Grad[ib][jb][k]+=grad*F_k[ib][jb][k];
			} else{
				tmp=Iq_exp[h0]-InC[i][j][k].real();
				energy+=tmp*tmp*Par*wei;
				grad=-2.0*Par*wei*tmp*scalePlot;
				Grad[i][j][k]+=2.0*grad*F_k[i][j][k];
			}
		}
	}

	F_k=Grad;

	return energy;
}
void Funktionell::Write(){
	std::string fileout="out_min.dat";
	std::ofstream fout(fileout.c_str(),ios::out);
	fout << std::scientific;
	for(auto it=Iq_c.begin();it!=Iq_c.end();it++){
		auto h0=it->first;
		fout << (h0+0.5)*dq<< " " << Iq_c[h0]*myScale<< " " << Iq_exp[h0]<<endl;
	}
}
double Funktionell::Energy2(double & scalePlot,double & Gscale,array3<Complex> & F_k){
	array3<Complex> Grad(nx,ny,nzp,sizeof(Complex));

	double energy{0};
	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
	,nfz{(nz % 2 == 0)? nz/2: nz/2+1};
	auto InC=Modulus(F_k);

	energy=0.0;
	Complex vt0,vt1;
	Grad=Complex{0,0};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		const size_t h0=it->first;
		vector<vector<int>> & vIdx=it->second;
		for(size_t o{0};o<vIdx.size();o++){
			int i=vIdx[o][XX];
			int j=vIdx[o][YY];
			int k=vIdx[o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				Grad[i][j][k]+=F_k[i][j][k];
				Grad[ib][jb][k]+=F_k[ib][jb][k];
				vt0=InC[i][j][k]+InC[ib][jb][k];
			} else{
				Grad[i][j][k]+=2.0*F_k[i][j][k];
				vt0=InC[i][j][k];
			}
			energy+=vt0.real();
		}
	}
//	double mw1,mw2,mw3,mw;
//	for(auto i=0;i<nx;i++){
//		int ia=(i<nfx)?i : i-nx;
//		size_t ib=i==0?0:nx-i;
//		for(auto j=0;j<ny;j++){
//			int ja=(j<nfy)?j : j-ny;
//			size_t jb=j==0?0:ny-j;
//			for(auto k=0;k<nzp;k++){
//				int ka=(k<nfz)?k : k-nz;
//				mw1=OC[XX][XX]*ia+OC[XX][YY]*ja+OC[XX][ZZ]*ka;
//				mw2=OC[YY][XX]*ia+OC[YY][YY]*ja+OC[YY][ZZ]*ka;
//				mw3=OC[ZZ][XX]*ia+OC[ZZ][YY]*ja+OC[ZZ][ZZ]*ka;
//				mw1=2.0*M_PI*mw1;
//				mw2=2.0*M_PI*mw2;
//				mw3=2.0*M_PI*mw3;
//				mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
//				int h0=static_cast<int>(mw/dq);
//				if(mw<qcut){
//					if(k != 0 && k != nzp-1){
//						Grad[i][j][k]+=F_k[i][j][k];
//						Grad[ib][jb][k]+=F_k[ib][jb][k];
//						vt0=InC[i][j][k]+InC[ib][jb][k];
//					} else{
//						Grad[i][j][k]+=2.0*F_k[i][j][k];
//						vt0=InC[i][j][k];
//					}
//					energy+=vt0.real();
//				}
//			}
//		}
//	}
	F_k=Complex{0,0};
	F_k=Grad;

	return energy;
}

double Funktionell::Deviate2(array3<Complex> & F_k){
	auto InC=Modulus(F_k);
	double energy{0};
	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};
	double mw1,mw2,mw3,mw;
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		const size_t h0=it->first;
		vector<vector<int>> & vIdx=it->second;
		vector<Dvect> & vV=mapVec[h0];
		for(size_t o{0};o<vIdx.size();o++){
			int i=vIdx[o][XX];
			int j=vIdx[o][YY];
			int k=vIdx[o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			mw1=vV[o][XX];
			mw2=vV[o][YY];
			mw3=vV[o][ZZ];
			mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
			if(mw<qcut){
				Complex v0=InC[i][j][k];
				if(k != 0 && k != nzp-1){
					Complex vt1=InC[ib][jb][k];
					v0=0.5*(v0+vt1);
				}
				int h0=static_cast<int>(mw/dq);
			}
		}
	}
	return energy;
}

array3<Complex> Funktionell::Modulus(array3<Complex> & ro_k){
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
Funktionell::~Funktionell() {
}

} /* namespace Funkll */

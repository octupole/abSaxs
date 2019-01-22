/*
 * Funktionell.cpp
 *
 *  Created on: Dec 13, 2018
 *      Author: marchi
 */

#include "Funktionell.h"
const double INEN{100000.0};
namespace Funkll {
struct pickPar{
	pickPar(double &Par, double & A0, double A1,map<size_t,double> & Ie,map<size_t,double> & Ic){
		double energy{0};
		auto i0=Ie.begin();
		for(auto it{Ie.begin()};it!=Ie.end();it++){
			auto h0=it->first;
			double tmp=Ie[h0]-Ic[h0]*A0-A1;
			double wei=1.0;
			energy+=tmp*tmp;
		}
//		Par=INEN/energy;
	Par=1.0/Ie[Ie.begin()->first];
	Par=1.0;
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
	Rge=exp->getRg();
	Rge*=Rge;
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
double Funktionell::EnergyR(double & myRg,array3<double> & F_r,array3<double> & gradR){
	double energy{0};
	/* Assume Cell axis are orthogonal */
	double PotR{0.2};
	double dx=co[XX][XX]/(double) Nx;
	double dy=co[YY][YY]/(double) Ny;
	double dz=co[ZZ][ZZ]/(double) Nz;
	Dvect c0;
	double norm{0};
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				double x0=o*dx;
				double y0=p*dy;
				double z0=q*dz;
				double density=F_r[o][p][q];
				c0[XX]+=x0*density;
				c0[YY]+=y0*density;
				c0[ZZ]+=z0*density;
				norm+=density;
			}
	c0/=norm;

	double Rg{0};
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				double x0=o*dx;
				double y0=p*dy;
				double z0=q*dz;
				double x1=c0[XX]-x0;
				double y1=c0[YY]-y0;
				double z1=c0[ZZ]-z0;
				double density=F_r[o][p][q];
				Rg+=(x1*x1+y1*y1+z1*z1)*density;
			}
	Rg/=norm;
	energy=PotR*(Rge-Rg)*(Rge-Rg);
	myRg=Rg<0?sqrt(-Rg):sqrt(Rg);
	double grad{-2.0*PotR*(Rge-Rg)};
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				double x0=o*dx;
				double y0=p*dy;
				double z0=q*dz;

				double x1=c0[XX]-x0;
				double y1=c0[YY]-y0;
				double z1=c0[ZZ]-z0;
				double one=x1*x1+y1*y1+z1*z1-Rg;
				one/=norm;
				gradR[o][p][q]+=one*grad;
			}

	return energy;
}
void Funktionell::ComputeIqc(array3<Complex> & F_k){
	Iq_c.clear();
	auto InC=Modulus(F_k);
	Complex vt0{0,0};
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
				vt0=0.5*InC[i][j][k]+0.5*InC[ib][jb][k];
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());
	}
}
double Funktionell::EnergyQ(double & AA_0,double & AA_1,array3<Complex> & F_k){
	A_1=0;
	A_0=AA_0;
//	A_0=AA_0;
//	A_1=AA_1;
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
				vt0=0.5*InC[i][j][k]+0.5*InC[ib][jb][k];
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());
	}
	static pickPar Once(Par,A_0,A_1,Iq_exp,Iq_c);
	Grad=Complex{0,0};
	energy=0;
	double scale{1.0/(double) mapIdx.size()};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		auto iqe=Iq_exp[h0];
		auto tmp=(Iq_exp[h0]-Iq_c[h0]*A_0-A_1);
		double wei=scale/iqe;
		energy+=tmp*tmp*Par*wei;

		double grad=-2.0*Par*wei*tmp*A_0/(double) mapIdx[h0].size();
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				Grad[i][j][k]+=grad*F_k[i][j][k]*0.5;
				Grad[ib][jb][k]+=grad*F_k[ib][jb][k]*0.5;
			} else{
				Grad[i][j][k]+=2.0*grad*F_k[i][j][k];
			}
		}
	}
	Grad[0][0][0]=Complex{0,0};
	F_k=Grad;
	return energy;
}
void Funktionell::scaleIq(array3<Complex> & F_k){
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
				vt0=0.5*(InC[i][j][k]+InC[ib][jb][k]);
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());
		auto scale=pow(Iq_exp[h0]/Iq_c[h0],1.0/4);
		for(size_t o{0};o<vIdx.size();o++){
			int i=vIdx[o][XX];
			int j=vIdx[o][YY];
			int k=vIdx[o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				F_k[i][j][k]*=Complex{scale,0};
				F_k[ib][jb][k]*=Complex{scale,0};
			} else{
				F_k[i][j][k]*=Complex{scale*scale,0};
			}
		}
	}
	Iq_c.clear();
	InC=Modulus(F_k);
	vt0=Complex{0,0};
	energy=0;
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
				vt0=0.5*(InC[i][j][k]+InC[ib][jb][k]);
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());
	}
}
double Funktionell::EnergyQ(double & AA_0,double & AA_1,array3<Complex> & F_k
		,array3<Complex> &DGrad){
	A_0=AA_0;
	A_1=AA_1;
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
				vt0=0.5*InC[i][j][k]+0.5*InC[ib][jb][k];
			} else{
				vt0=InC[i][j][k];
			}
			Iq_c[h0]+=vt0.real();
		}
		Iq_c[h0]/=static_cast<double>(vIdx.size());

	}
	static pickPar Once(Par,A_0,A_1,Iq_exp,Iq_c);
	Grad=Complex{0,0};
	DGrad=Complex{0,0};
	energy=0;
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		auto tmp=(Iq_exp[h0]-Iq_c[h0]*A_0-A_1);
		double wei=1.0;
		energy+=tmp*tmp*Par*wei;
		double dtmp=-2.0*tmp*A_0/(double) mapIdx[h0].size();
		double ddtmp=-2.0*Par*wei*A_0*A_0/(double) (mapIdx[h0].size()*mapIdx[h0].size());
		double grad=-2.0*Par*wei*tmp*A_0/(double) mapIdx[h0].size();
		//		dA_0+=-2.0*Par*wei*tmp*Iq_c[h0];
		//		dA_1+=-2.0*Par*wei*tmp;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			size_t ib=i==0?0:nx-i;
			size_t jb=j==0?0:ny-j;
			if(k != 0 && k != nzp-1){
				auto f_k=F_k[i][j][k]*0.5;
				auto f_kb=F_k[ib][jb][k]*0.5;
				Grad[i][j][k]+=grad*f_k;
				Grad[ib][jb][k]+=grad*f_kb;
			} else{
				auto f_k=F_k[i][j][k];
				Grad[i][j][k]+=2.0*grad*f_k;
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
		fout << (h0+0.5)*dq<< " " << Iq_c[h0]*A_0+A_1<< " " << Iq_exp[h0]<<endl;
	}
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

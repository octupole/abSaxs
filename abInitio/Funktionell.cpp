/*
 * Funktionell.cpp
 *
 *  Created on: Dec 13, 2018
 *      Author: marchi
 */

#include "Funktionell.h"

namespace Funkll {

Funktionell::Funktionell(RhoSaxs * ro_out, RhoSaxs * ro_in,SaxsData * exp){
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
	for(auto it{qft.begin()};it!=qft.end();it++){
		const int h0=it->first;
		Iq_exp[h0]=qft[h0]/static_cast<double>(nqft[h0]);
		qcut=h0*dq+0.5*dq;
	}

	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};
	double mw1,mw2,mw3,mw;
	size_t MM{0};
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
				if(mw<qcut){
					vInt[h0].push_back(vector<int>{i,j,k});
					vDble[h0].push_back(Dvect{mw1,mw2,mw3});
					MM++;
				}
			}
		}
	}
}

double Funktionell::Deviate(array3<Complex> & F_k){
	auto InC=Modulus(F_k);
	double energy{0};
	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};
	double mw1,mw2,mw3,mw;
	for(auto it=vInt.begin();it != vInt.end();it++){
		const size_t h0=it->first;
		vector<vector<int>> & vIdx=it->second;
		vector<Dvect> & vV=vDble[h0];
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
			Complex v0=InC[i][j][k];
			if(k != 0 && k != nzp-1){
				Complex vt1=InC[ib][jb][k];
				v0=0.5*(v0+vt1);
			}
			Iq_c[h0].first+=v0.real();
			Iq_c[h0].second++;
		}
	}
	auto it_e=Iq_exp.begin();
	auto it_c=Iq_c.begin();


	for(auto it=Iq_c.begin();it!=Iq_c.end();it++){
		auto h0=it->first;
		Iq_c[h0].first/=static_cast<double>(Iq_c[h0].second);
	}
	this->Scaling=it_e->second/it_c->second.first;

	for(auto it=Iq_exp.begin();it!=Iq_exp.end();it++){
		auto h0=it->first;
		auto tmp=Iq_exp[h0]-Iq_c[h0].first*Scaling;
		energy+=tmp*tmp;
		cout << h0*dq+0.5*dq << " "  <<std::scientific << Scaling*Iq_c[h0].first<< " " <<Iq_exp[h0]<<endl;
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

/*
 * ABSaxs.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: marchi
 */

#include "ABSaxs.h"
namespace abinit{
class LBFGSWrapper{
	static ABSaxs * inst;
public:
	static void setInstance(ABSaxs * x){inst=x;}
	static void myFunc(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr){
		return inst->minLBFGS(x,func,grad,ptr);
	}

};
class LBFGSWrapper_C{
	static ABSaxs * inst;
public:
	static void setInstance(ABSaxs * x){inst=x;}
	static void myFunc(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr){
		return inst->minLBFGS_C(x,func,grad,ptr);
	}

};
ABSaxs * LBFGSWrapper::inst=nullptr;
ABSaxs * LBFGSWrapper_C::inst=nullptr;
void ABSaxs::testDensity(array3<double> & F){
	double avg{0},mMax{-1}, mMin{1};
	int M{0};
	double intgr{0.0};
	for(size_t o{0};o<F.Nx()*F.Ny()*F.Nz();o++){
		double tmp=*(&F[0][0][0]+o);
		if(fabs(tmp) <1.0e-5) continue;
		avg+=tmp;
		if(mMax <tmp)mMax=tmp;
		if(mMin >tmp)mMin=tmp;
		M++;
		}
	cout <<" Total electrons "<< avg*dvol << " Max= " << mMax << " Min = " <<mMin<< " " << M<<endl;;
}
ABSaxs::ABSaxs(uint nx, uint ny, uint nz, double SupCell): grid_b{nx,ny,nz},
		nx{nx},ny{ny},nz{nz},nzp{nz/2+1},SuperCell{SupCell} {
		}
void ABSaxs::setUpRho(SaxsData * exp){
	double Rd=exp->Rd();
	Rho_in->initDensity(Rd);
	Rho_s->PartialCopy(*Rho_in);
	myFuncx=new Funkll::Funktionell(Rho_s,Rho_in,exp);
}
void ABSaxs::setUpRho(SaxsData * exp,array3<double> & F_r){
	try{
	if(!bCellCalled) throw string("Cannot call setUpRho before setUpCell.");
	}catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	F_r/=dvol;
	for(size_t o{0};o<F_r.Nx();o++)
		for(size_t p{0};p<F_r.Ny();p++)
			for(size_t q{0};q<F_r.Nz();q++)
				(*Rho_s)[0][o][p][q]=F_r[o][p][q];
	Rho_in->PartialCopy(*Rho_s);
	myFuncx=new Funkll::Funktionell(Rho_s,Rho_in,exp);
}

void ABSaxs::setUpCell(SaxsData * exp){
	this->bCellCalled=true;
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

	Rho_in=new abInitioRho::RhoSaxs(grid_a[XX],grid_a[YY],grid_a[ZZ],co);
	Rho_s=new abInitioRho::RhoSaxs(grid_b[XX],grid_b[YY],grid_b[ZZ],CO);

	ptrdiff_t alloc_local, local_n0, local_0_start;

     alloc_local = fftw_mpi_local_size_3d(nx, ny, nzp, MPI_COMM_WORLD,
                                            &local_n0, &local_0_start);

  	F_k.Allocate(grid_b[XX],grid_b[YY],nzp,Calign);
 	F_k0.Allocate(grid_b[XX],grid_b[YY],nzp,Calign);
 	DGrad.Allocate(grid_b[XX],grid_b[YY],nzp,Calign);
	F_r.Allocate(grid_b[XX],grid_b[YY],nzp*2,Ralign);
	GradR.Allocate(grid_b[XX],grid_b[YY],nzp*2,Ralign);
	I_k.Allocate(grid_b[XX],grid_b[YY],nzp,Calign);
	cout << "Run with Rho_inner ="<< Nx<<" "<< Ny<<" "<< Nz<<" Resolution "<< co[XX][XX]/(float) Nx<<endl;
	cout << co;
	cout << "Run with Rho_outer ="<< nx<<" "<< ny<<" "<< nz<<" "<<endl;
	cout << CO;
}
void ABSaxs::Run(){
	Minimize();
//	this->testGradient();
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
void ABSaxs::minLBFGS(const real_1d_array &x, double & energy, real_1d_array &grad, void *ptr){
	size_t M{0};
	static int NN=0;
	F_r=0.0;
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				F_r[o][p][q]=x[M++];
				GradR[o][p][q]=0;
			}
	Forward3->fft(F_r,F_k);
	energy=0;
	energy+=myFuncx->EnergyQ(A_0,A_1,F_k);
	Backward3->fft(F_k,GradR);
	double myRg{0};
	energy+=myFuncx->EnergyR(myRg,F_r,GradR);
	M=0;
	double msd{-2};
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				grad[M]=GradR[o][p][q];
				if(msd < fabs(grad[M])) msd=fabs(grad[M]);
				M++;
			}
	cout << std::scientific;
	cout << "Step "<< NN++<< " Energy = "<<energy << " Grad = " << fabs(msd)
			<<" Rg = "<< std::fixed<<myRg<<endl;
}
void ABSaxs::minLBFGS_C(const real_1d_array &x, double & energy, real_1d_array &grad, void *ptr){
	size_t M{0};
	static int NN=0;
	F_k0=Complex{0.0,0.0};

	auto mapIdx=myFuncx->getIdx();
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			F_k0[i][j][k]=Complex{x[M++],x[M++]};
		}
	}
	energy=myFuncx->EnergyQ(A_0,A_1,F_k0);
	M=0;
	double ddd{1};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			grad[M]=2.0*F_k0[i][j][k].real();
			grad[M+1]=2.0*F_k0[i][j][k].imag();
			M+=2;
		}
	}

	double msd{-2};
	for(size_t o{0};o<M;o++){
		if(msd< fabs(grad[o])) msd=fabs(grad[o]);
	}
	cout << std::scientific;

	cout << "Step "<< NN++<< " Energy = "<<energy << " Grad = " << fabs(msd)<<endl;
}
void ABSaxs::subtract(array3<double> & F){
	for(size_t o{0};o<F.Size();o++){
		auto & tmp=*(&F[0][0][0]+o);
		if(tmp < wDensity){
			tmp=0.0;
		}else{
			tmp-=wDensity;
		}
	}

}
void ABSaxs::Minimize(){
    double epsg = 0.01;
    double epsf = 0;
    double epsx = 0;
    ae_int_t maxits = 0;

    minlbfgsreport rep;
	real_1d_array x0,grad;
	x0.setlength(Nx*Ny*Nz);
	grad.setlength(Nx*Ny*Nz);
	/*
	 * Arrays are destroyed while performing plan!!!
	 */

	Forward3 =new Pfftwpp::Prcfft3d(nx,ny,nz,F_r,F_k);
	Backward3=new Pfftwpp::Pcrfft3d(nx,ny,nz,F_k,GradR);

	Rho_s->copyOut(F_r);
	this->subtract(F_r);
	this->_radius();
	F_r*=dvol;
	Forward3->fft(F_r,F_k);
	myFuncx->ComputeIqc(F_k);
	auto Iqc=myFuncx->getIqc();
	auto Iqe=myFuncx->getIqe();
	dq=myFuncx->getMydq();

	this->fitSaxs(Iqe,Iqc);
	A_0=1;
	A_1=0;

	Rho_s->copyOut(F_r);
	this->subtract(F_r);
	F_r*=dvol;
	size_t M{0};
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				x0[M++]=F_r[o][p][q];
			}
	LBFGSWrapper::setInstance(this);

    minlbfgscreate(5, x0, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    minlbfgssetxrep(state, true);

    alglib::minlbfgsoptimize(state, &LBFGSWrapper::myFunc);
    minlbfgsresults(state, x0, rep);

	(*Rho_s)=0;
	this->_radius();
    Rho_s->copyIn(F_r);
    Rho_s->WriteIt();
    myFuncx->Write();
}
void ABSaxs::Minimize_C(){
    double epsg = 0.001;
    double epsf = 0;
    double epsx = 0;
    ae_int_t maxits = 20;

    minlbfgsreport rep;
	real_1d_array x0,grad;
	/*
	 * Arrays are destroyed while performing plan!!!
	 */

	Forward3 =new Pfftwpp::Prcfft3d(nx,ny,nz,F_r,F_k);
	Backward3=new Pfftwpp::Pcrfft3d(nx,ny,nz,F_k,F_r);


	Rho_s->copyOut(F_r);
	this->_radius();
	this->subtract(F_r);
//	this->testDensity(F_r);exit(1);

	F_r*=dvol;
	Forward3->fft(F_r,F_k);
	myFuncx->ComputeIqc(F_k);
	auto Iqc=myFuncx->getIqc();
	auto Iqe=myFuncx->getIqe();
	auto mapIdx=myFuncx->getIdx();
	dq=myFuncx->getMydq();

	this->fitSaxs(Iqe,Iqc);
	size_t M{0};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			M+=2;
		}
	}

	x0.setlength(M);
	grad.setlength(M);
	dgrad.setlength(M);
//	alglib::minlbfgssetprecdiag(state,dgrad);


	M=0;
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			x0[M++]=F_k[i][j][k].real();
			x0[M++]=F_k[i][j][k].imag();
		}
	}
	LBFGSWrapper_C::setInstance(this);

    minlbfgscreate(6, x0, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    minlbfgssetxrep(state, true);

    alglib::minlbfgsoptimize(state, &LBFGSWrapper_C::myFunc);
    minlbfgsresults(state, x0, rep);
	M=0;
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			F_k[i][j][k]=Complex{x0[M++],x0[M++]};
		}
	}
	Backward3->fftNormalized(F_k,F_r);


	Rho_s->copyIn(F_r);
    Rho_s->WriteIt();
	this->_radius();
    myFuncx->Write();
}
void ABSaxs::_radius(){
	Dvect c0;
	double norm{0};
	size_t N{0};
//	for(size_t o{0};o<nx;o++)
//		for(size_t p{0};p<ny;p++)
//			for(size_t q{0};q<nzp*2;q++){
//				if(fabs(F_r[o][p][q])< 0.09)F_r[o][p][q]=0;
//			}
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
	for(size_t o{0};o<3;o++)
		c0[o]/=norm;
	double Rg{0};
	N=0;
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
				N++;

			}
	Rg/=norm;
	cout << "Rg = "<<sqrt(Rg) <<" " << norm<<endl;

}
void ABSaxs::testMinim(){

}
void ABSaxs::testGradient(){
	real_1d_array x0,grad;
	x0.setlength(Nx*Ny*Nz);
	grad.setlength(Nx*Ny*Nz);
	/*
	 * Arrays are destroyed while performing plan!!!
	 */

	Forward3 =new Pfftwpp::Prcfft3d(nx,ny,nz,F_r,F_k);
	Backward3=new Pfftwpp::Pcrfft3d(nx,ny,nz,F_k,GradR);
	Rho_s->copyOut(F_r);
	F_r*=dvol;
	Forward3->fft(F_r,F_k);
	myFuncx->ComputeIqc(F_k);
	auto Iqc=myFuncx->getIqc();
	auto Iqe=myFuncx->getIqe();
	dq=myFuncx->getMydq();

	this->fitSaxs(Iqe,Iqc);
	Rho_s->copyOut(F_r);
	F_r*=dvol;

	size_t M{0};
	for(size_t o{0};o<Nx;o++)
		for(size_t p{0};p<Ny;p++)
			for(size_t q{0};q<Nz;q++){
				x0[M++]=F_r[o][p][q];
			}


	double c=F_r[2][2][3];
	double delta=0.002,eps=0.002;
	double Beg{-delta*0.5};
	double End{delta*0.5};
	vector<double> ee{0,0,0};
	M=0;
	void * myptr{nullptr};
	auto time1=MPI_Wtime();
	for(auto o{Beg};o<=End;o+=delta*0.5){
		F_r[2][2][3]=c+o;
		x0[2*Nz*Ny+2*Nz+3]=c+o;
		double E;
		this->minLBFGS(x0, E,grad, myptr);
		ee[M++]=E;
		cout << " E =" << E;
		cout <<" Grad " << GradR[2][2][3]<<endl;
	}
	cout << (ee[2]-ee[0])/delta<< endl;
	auto time2=MPI_Wtime();

	cout << time2-time1<<endl;
}
void ABSaxs::testGradient_C(){
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution (0.0,1.0);

	real_1d_array x0,grad;
	Forward3 =new Pfftwpp::Prcfft3d(nx,ny,nz,F_r,F_k);
	Backward3=new Pfftwpp::Pcrfft3d(nx,ny,nz,F_k,F_r);
	array3<Complex> G(F_k);
	Rho_s->copyOut(F_r);
	for(auto o=0;o<F_r.Size();o++)
		*(&F_r[0][0][0]+o)=distribution(generator);
	F_r*=dvol;
	Forward3->fft(F_r,F_k);
	myFuncx->ComputeIqc(F_k);
	auto Iqc=myFuncx->getIqc();
	auto Iqe=myFuncx->getIqe();
	auto mapIdx=myFuncx->getIdx();
	dq=myFuncx->getMydq();

	this->fitSaxs(Iqe,Iqc);

	cout << F_k[2][2][3].real() <<endl;

	double c1=F_k[2][2][3].imag();
	double c=F_k[2][2][3].real();

	double delta=0.2,eps=0.002;
	double Beg{-delta*0.5};
	double End{delta*0.5};
	vector<double> ee{0,0,0};
	int M=0;
	void * myptr{nullptr};

	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			M+=2;
		}
	}
	x0.setlength(M);
	grad.setlength(M);
	auto time1=MPI_Wtime();
	int MM{0};
	M=0;
	int key0{0},key1{0};
	for(auto it=mapIdx.begin();it != mapIdx.end();it++){
		auto h0=it->first;
		for(size_t o{0}; o<mapIdx[h0].size();o++){
			size_t i=mapIdx[h0][o][XX];
			size_t j=mapIdx[h0][o][YY];
			size_t k=mapIdx[h0][o][ZZ];
			if(i == 2 && j == 2 && k == 3) {
				cout << M+1 << endl;
			}
			x0[M++]=F_k[i][j][k].real();
			x0[M++]=F_k[i][j][k].imag();
			if(i == 0 && j == 2 && k == 12) {
				key0=M-2;key1=M-1;
			}
		}
	}
		cout << key0<<endl;
		cout << x0[key1]<< " " << F_k[0][2][12].imag()<<endl;
	for(auto o{Beg};o<=End;o+=delta*0.5){
		F_k[0][2][12]=Complex{c+o,c1};
		x0[key1]=c+o;
		cout<< "o = " << o ;
		double E;
		this->minLBFGS_C(x0, E,grad, myptr);
		ee[MM++]=E;
		cout << " E =" << E;
		cout <<" Grad " << grad[key1] <<endl;
	}
	cout << (ee[2]-ee[0])/delta<< endl;
	auto time2=MPI_Wtime();

	cout << time2-time1<<endl;
//	for(size_t o{0};o<Nx;o++)
//		for(size_t p{0};p<Ny;p++)
//			for(size_t q{0};q<Nz;q++){
//				if(GradR[o][p][q] > 1.0e-2)
//					cout << o<< " " <<p<< " " <<q<< " " <<GradR[o][p][q]<<endl;
//			}


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
void ABSaxs::fitSaxs(map<size_t,double> & Iqe, map<size_t,double> & Iqc){
	real_2d_array xy;
	real_1d_array S;
	xy.setlength(Iqe.size(),2);
	S.setlength(Iqe.size());
	int o{0};
	for(auto op{Iqe.begin()};op!=Iqe.end();op++){
		auto h0=op->first;
		xy[o][0]=Iqc[h0];
		xy[o][1]=Iqe[h0];
		S[o]=xy[o][1];
		o++;
	}
	ae_int_t info;
    ae_int_t nvars;
    linearmodel model;
    lrreport rep;
    real_1d_array c;
    lrbuilds(xy,S,Iqe.size(), 1, info, model, rep);

    lrunpack(model, c, nvars);

    A_0=c[0];
    A_1=c[1];
//    o=0;
//    for(auto op{Iqe.begin()};op!=Iqe.end();op++){
//    	auto h0=op->first;
//    	double y=(c[0]*Iqc[h0]+c[1]);
//    	cout << (h0+0.5)*dq <<" "<< y<< " " << Iqc[h0]<< " "<<Iqe[h0]<<endl;
//    	o++;
//    }
// 	exit(1);
}

ABSaxs::~ABSaxs() {
	// TODO Auto-generated destructor stub
}
}

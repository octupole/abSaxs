/*
 * NewMPI.h
 *
 *  Created on: Jan 28, 2016
 *      Author: marchi
 */

#ifndef LIBTRAJ_NEWMPI_H_
#define LIBTRAJ_NEWMPI_H_

#include <iostream>
#include <vector>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <complex>


namespace MPI{
template <bool MMPI>
class Intracom;
}
#ifdef HAVE_MPI
#include <mpi.h>
#define OMPI_SKIP_MPICXX
const bool mMPI{true};
#else
const bool mMPI{false};

static void * MPI_STATUS_IGNORE;
static void * MPI_COMM_WORLD;
static void * MPI_INT;
static void * MPI_FLOAT;
static void * MPI_SUM;
static void * MPI_DOUBLE;
static void * MPI_DOUBLE_COMPLEX;
static void * MPI_IN_PLACE;
static void * MPI_CHAR;
using MPI_Comm = void *;
using MPI_Datatype = void *;
using MPI_Op = void *;
inline static void MPI_Init(...){};
inline static void MPI_Finalize(){};
inline static void MPI_Comm_rank(...){};
inline static void MPI_Comm_size(...){};
inline static void MPI_Barrier(...){};
inline static void MPI_Bcast(...){};
inline static void MPI_Reduce(...){};
inline static void MPI_Gather(...){};
inline static void MPI_Allgather(...){};
inline static void MPI_Send(...){};
inline static void MPI_Recv(...){};

inline static double MPI_Wtime(){
	timeval tim;
	double temp=tim.tv_sec+(tim.tv_usec/1000000.0);
	return temp;
	}

#endif


using std::vector;
using std::string;
using std::cout;
using std::endl;

namespace MPI{
const size_t zero{0};
const size_t one{1};

template <>
class Intracom<false>{
public:
	Intracom(void *){};
	int Get_size(){return one;}
	int Get_rank(){return zero;}
	void Barrier(){};
	void Reduce(...){};
	void Bcast(...){};
	void Gather(...){};
	void Allgather(...){};
	void Send(...){};
	void Recv(...){};
	void set(...){};
	Intracom(){};

	virtual ~Intracom(){};
};
template <>
class Intracom<true>{
	MPI_Comm comm;
public:
	Intracom(MPI_Comm x):comm{x}{};
	void set(MPI_Comm x){
		comm=x;
	};
	int Get_rank(){
		int rank{zero};
		try{
			if(MPI_Comm_rank(comm,&rank) != 0) throw string("Error obtaining rank!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
		return rank;
		}
	int Get_size(){
		int size{one};
		try{
			if(MPI_Comm_size(comm,&size) != 0) throw string("Error obtaining size!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
		return size;
	}

	void Barrier(){
		MPI_Barrier(comm);
	};
	void Reduce(void * sendbuf, void * recvbuf, int count, MPI_Datatype datatype,MPI_Op op,int root){
		try{
			int ierr=MPI_Reduce(sendbuf,recvbuf,count,datatype,op,root,comm);
			if(ierr) throw string("Error from MPI_Reduce!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
		
	};
	void Bcast(void* buffer, int count, MPI_Datatype datatype, int root){
		try{
			int ierr=MPI_Bcast(buffer,count, datatype, root,comm);
			if(ierr) throw string("Error from MPI_Bcast!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
	};
	void Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,void* recvbuf, int recvcount, MPI_Datatype recvtype, int root){
		try{
			int ierr=MPI_Gather(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm);
			if(ierr) throw string("Error from MPI_Gather!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}

	};
	void Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,void* recvbuf, int recvcount, MPI_Datatype recvtype){
		try{
			int ierr=MPI_Allgather(sendbuf, sendcount, sendtype,recvbuf, recvcount,recvtype,comm);
			if(ierr) throw string("Error from MPI_Allgather!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
	};
	void Send(void* buf, int count, MPI_Datatype datatype, int dest,int tag){
		try{
			int ierr=MPI_Send(buf, count, datatype,dest,tag,comm);
			if(ierr) throw string("Error from MPI_Send!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
	};
	void Recv(void* buf, int count, MPI_Datatype datatype, int source,int tag){
		try{
			int ierr=MPI_Recv(buf, count,datatype,source,tag,comm,MPI_STATUS_IGNORE);
			if(ierr) throw string("Error from MPI_Recv!");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
	};
	virtual ~Intracom(){}
};
}




namespace Parallel{
typedef std::complex<double> Complex;
const bool DEBUG{false};
class NewMPI {
	bool is_parallel{mMPI};
	MPI::Intracom<mMPI> comm{MPI_COMM_WORLD};
public:
	NewMPI(int & argc,char ** & argv){
		MPI_Init(&argc,&argv);
		if(is_parallel){
			if(comm.Get_size()){
				if(comm.Get_rank()){
					std::cout.setstate(std::ios_base::badbit);
					fclose(stdout);
					fclose(stderr);
				}
				std::cout << "\n" << std::endl;
				std::cout << " ------ Parallel run with " << comm.Get_size() << " CPUS " << std::endl;
			}
		}
	};
	NewMPI(){
		MPI_Init(nullptr,nullptr);
		if(is_parallel){
			if(comm.Get_size()){

				if(comm.Get_rank()){
					std::cout.setstate(std::ios_base::badbit);
					fclose(stdout);
					fclose(stderr);
				}
				std::cout << "\n" << std::endl;
				std::cout << " ------ Parallel run with " << comm.Get_size() << " CPUS " << std::endl;
				is_parallel=true;
			}
		}
	};
	size_t Get_Rank(){
		return comm.Get_rank();
	}
	size_t Get_Size(){
		return comm.Get_size();
	}
	bool AmI_Parallel(){return is_parallel;};

	template <class T>
	void Broadcast(T *, const int);

	template <class T>
	void ReduceSum(T *, const int);

	template <class T>
	void ReduceSum(vector<T> &);

	template <class T>
	void Gather(vector<T> & , vector<T> &);

	template <typename T>
	void AllGather(int, T *, T * );

	void Barrier(){comm.Barrier();}
	double Time(){return MPI_Wtime();};

	template <typename T>
	void Send(int dest, int dim, int tag, T * sbuffer);

	template <typename T>
	void Recv(int source,int dim, int tag, T * rbuffer);

	MPI_Comm * Communicator(){ return &comm;}
	void Finalize(){MPI_Finalize();};
	virtual ~NewMPI(){MPI_Finalize();};
};


template<>
inline void NewMPI::Gather<double>(vector<double> & bufferIn, vector<double> & bufferOut){
	int nIn=static_cast<int> (bufferIn.size());
	int nbufferIn,nbufferOut;
	nbufferIn=nIn;
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(&nIn,NULL,1,MPI_INT,MPI_SUM,0);
	else comm.Reduce(MPI_IN_PLACE,&nIn,1,MPI_INT,MPI_SUM,0);
	comm.Bcast(&nIn,1,MPI_INT,0);
	nbufferOut=nIn;
	bufferOut.clear();
	bufferOut=vector<double>(nbufferOut,0.0);
	comm.Gather(&bufferIn[0],nbufferIn,MPI_DOUBLE,&bufferOut[0],nbufferIn,MPI_DOUBLE,0);
}

template<>
inline void NewMPI::Broadcast<double>(double * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI_DOUBLE,0);
};

template<>
inline void NewMPI::Broadcast<int>(int * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI_INT,0);
	};
template<>
inline void NewMPI::Broadcast<float>(float * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI_FLOAT,0);
};

template<>
inline void NewMPI::Broadcast<char>(char * buffer,const int nbuffer){
	comm.Bcast(buffer,nbuffer,MPI_CHAR,0);
	};

template<>
inline void NewMPI::ReduceSum<float>(float * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI_FLOAT,MPI_SUM,0);
	else comm.Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_FLOAT,MPI_SUM,0);
	comm.Barrier();
	};
template<>
inline void NewMPI::ReduceSum<double>(double * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI_DOUBLE,MPI_SUM,0);
	else comm.Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_DOUBLE,MPI_SUM,0);
	comm.Barrier();
};
template<>
inline void NewMPI::ReduceSum<Complex>(Complex * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI_DOUBLE_COMPLEX,MPI_SUM,0);
	else comm.Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_DOUBLE_COMPLEX,MPI_SUM,0);
	comm.Barrier();
};
template<>
inline void NewMPI::ReduceSum<double>(vector<double> & buffer0){
	const int nbuffer=buffer0.size();
	double * buffer;
	buffer=new double[nbuffer];
	for(int n=0;n<nbuffer;n++)
		buffer[n]=buffer0[n];
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI_DOUBLE,MPI_SUM,0);
	else comm.Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_DOUBLE,MPI_SUM,0);
	for(int n=0;n<nbuffer;n++)
		buffer0[n]=buffer[n];
	comm.Barrier();
	};
template<>
inline void NewMPI::ReduceSum<int>(int * buffer,const int nbuffer){
	comm.Barrier();
	if(comm.Get_rank()) comm.Reduce(buffer,NULL,nbuffer,MPI_INT,MPI_SUM,0);
	else comm.Reduce(MPI_IN_PLACE,buffer,nbuffer,MPI_INT,MPI_SUM,0);
	};
template<>
inline void NewMPI::AllGather<int>(int n, int * sendbuff, int * recvbuf){
	comm.Allgather(sendbuff,n, MPI_INT, recvbuf, n, MPI_INT);
}
template<>
inline void NewMPI::AllGather<double>(int n, double * sendbuff, double * recvbuf){
	comm.Allgather(sendbuff,n, MPI_DOUBLE, recvbuf, n, MPI_DOUBLE);
}
template<>
inline void NewMPI::Send<float>(int dest, int nbuffer, int tag, float * x){
	comm.Send(x,nbuffer,MPI_FLOAT,dest,tag);
}
template<>
inline void NewMPI::Send<double>(int dest, int nbuffer, int tag, double * x){
	comm.Send(x,nbuffer,MPI_DOUBLE,dest,tag);
}
template<>
inline void NewMPI::Recv<float>(int source, int nbuffer, int tag, float * x){
	comm.Recv(x,nbuffer,MPI_FLOAT,source,tag);
}
template<>
inline void NewMPI::Recv<double>(int source, int nbuffer, int tag, double * x){
	comm.Recv(x,nbuffer,MPI_DOUBLE,source,tag);
}
} /* namespace Parallel */


#endif /* LIBTRAJ_NEWMPI_H_ */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <string>
#include <sstream>
#include "mpi.h"

using namespace std;


class vector
{
public:
	double *val;
	int length;
	vector():length(0){};
	vector(const int &t_n){
		length=t_n;
		val= new double[length];
	};	
	~vector(){
	length=0;
	delete[] val;
	}
	void print(){
		for(int i=0;i<length;i++){
			cout<<val[i]<<endl;
		}
	}

};
class matrix{
public:int nnz,nrows,ncols,*indi,*indj;
	   double *val;
		//Default Constructor
		matrix():nrows(0),ncols(0),nnz(0){};
		//Secondary Constructor
		//matrix(const int &t_nrows,const int &t_ncols,const int &t_nnz);
		~matrix(){delete[] val; delete[] indj; delete[] indi; ncols=0;nrows=0;nnz=0;};
		matrix(const int &t_nrows,const int &t_ncols,const int &t_nnz){
	
	nrows=t_nrows;
	ncols=t_ncols;
	nnz=t_nnz;
	val = new double[nnz];
	indj = new int[nnz];
	indi = new int[nrows+1];

	
	
}
	void print(){
		for(int i=0;i<nrows;i++){
            cout<<" ["<<indi[i]<<" , "<<indj[i]<<" ]";
			cout<<val[i]<<endl;	
		}
	}
};


void SpMV(vector &y,const matrix &M, const vector &x,const int &sumnnzm){
	double sum=0;
	#pragma omp parallel for private(i) lastprivate(sum)	
	for(int i = 0; i<M.nrows; i++){
		sum=0;	
			
		for (int j=(M.indi[i]-sumnnzm);j<(M.indi[i+1]-sumnnzm);j++){			
			sum+=M.val[j]*x.val[M.indj[j]];			
		}
		y.val[i]=sum;		
	}
}

void VVadd(vector &y,const vector &a, const vector &b){
	#pragma omp parallel for private(i)
	for(int i = 0 ; i<y.length; i++){
		y.val[i]=a.val[i]+b.val[i];
	}
}

void VVdic(vector &y,const vector &a, const vector &b){
	#pragma omp parallel for private(i) 
	for(int i = 0 ; i<y.length; i++){
		y.val[i]=a.val[i]-b.val[i];
	}
}

void VVmul(double &y,const vector &a, const vector &b){
        double sum=0;
	#pragma omp parallel for private(i) lastprivate(sum)
	for(int i = 0 ; i<a.length; i++){
		sum=sum+a.val[i]*b.val[i];
	}
	y=sum;
}

void CVmul(vector &y,double &a, const vector &b){
	#pragma omp parallel for private(i)
	for(int i = 0 ; i<y.length; i++){
		y.val[i]=a*b.val[i];
	}
}

void inverse(matrix &K,const matrix &A){
	#pragma omp parallel for private(i)
	for (int i=0; i<A.nrows; i++)
	{
		K.val[i]=1/A.val[i];
		K.indi[i]=A.indi[i];
		K.indj[i]=A.indj[i];	
	}
	K.indi[A.nrows]=A.indi[A.nrows];	
}

double norm(double &y,const vector &x){
	double sum=0;
	for (int i=0; i<x.length; i++){
		sum=sum+pow(x.val[i],2.0);
	}
	y=sqrt(sum);
}

void Jacobi(matrix &K,const matrix &A){ 
	for (int i=0; i<A.nrows; i++)
	{	K.indi[i]=i;	
		for(int j=A.indi[i]; j<A.indi[i+1]; j++){
			if (A.indj[j]==i){
				K.val[i]=A.val[j];
				K.indj[i]=A.indj[j];								
			}
		}
	}
	K.indi[A.nrows]=A.nrows;
}



int main(int argc, char **argv){

	int rank,procs,ncolsm,nrowsm,nozerosm,fnzm,sizeindim,fsearchm,lsearchm,sumnnzm;
	int *ncols, *nrows, *nozeros, *fnz, *lnz, *sizeindi, *fnzi, *fsearch, *lsearch,*sumnnz;
	
	//pername ton pinaka
	std::ifstream file("A.txt");
        std::string str; 
    	int k=0,nnz=0,size=0;
   	int n; // proti grammh - size
    	getline(file, str);    
    	std::stringstream stream(str);
    	stream>>n;
    	nnz=n;
    	//cout<<nnz<<endl;
    	getline(file, str);   // deuterh grammh gia to size  
    	std::stringstream stream2(str);
    	stream2>>n;
    	size=n;
    	//cout<<size<<endl;
    	matrix A(size,size,nnz); // dhlwnoume ton pinaka 
    	getline(file, str);   // trith grammh gia to val[]  
    	std::stringstream stream3(str);   
    	double vi;  
    	int i=0;      
    	while(stream3>>vi){                
         	if(!stream) break;     
         	A.val[i]=vi;
         	i++;                                          
    	}  
    	getline(file, str);   // 4h grammh gia to indj[]  
    	std::stringstream stream4(str);        
    	i=0;      
    	while(stream4>>n){                
         	if(!stream) break;     
         	A.indj[i]=n;
         	i++;                                          
    	} 
    	getline(file, str);   // 5h grammh gia to indi[]  
    	std::stringstream stream5(str);        
    	i=0;      
    	while(stream5>>n){                
         	if(!stream) break;     
         	A.indi[i]=n;
         	i++;                                          
    	}         
	std::ifstream file2("b.txt"); // dhmiourgoume to vector 
    	std::string str2; 
    	getline(file2, str2);
    	std::stringstream streamB1(str2);
    	streamB1>>n;
    	size=n;
    	vector b(size);
    	getline(file2, str2);   // fortonoume tis times sto vector mas   
    	std::stringstream streamB2(str2);        
    	i=0;      
    	while(streamB2>>vi){                
         	if(!stream) break;     
         	b.val[i]=vi;
         	i++;                                          
    	}
	

	

	//orizoume tis paralliles metavlites kai pinakes
	
	MPI_Status status;
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	double time = MPI_Wtime();
	
	
if (rank==0){
	ncols= new int [procs];
	nrows= new int [procs];
	nozeros=new int [procs];
	fnz= new int [procs];
	sizeindi= new int [procs];
	fnzi = new int [procs];
	fsearch= new int [procs];
	lsearch= new int [procs];
	
	sumnnz= new int [procs];
	int sum=0,sum2=0;

	for (int i=0; i<procs; i++){sumnnz[i]=0;}
	fnzi[0]=0;
	for (int i=0; i<procs; i++){
		ncols[i]=size;
		nrows[i]=size/procs;
		nozeros[i]=A.indi[(i+1)*nrows[i]]-A.indi[i*nrows[i]];
		sum+=nozeros[i];
		fnz[i]=sum-nozeros[i];
		
		
		if (i>0) sumnnz[i]=nozeros[i-1];	

		fsearch[i]=i*nrows[i];
		lsearch[i]=fsearch[i]+nrows[i]-1;
		
		sizeindi[i]=nrows[i]+1;
		fnzi[i+1]=nrows[i];
	}
	
	int count=0;
	int k=0;
	for (int i=0; i<(A.nrows); i++){
		count=count+1;
		if (A.indi[i]==fnz[k+1]){
			sizeindi[k]=count;
			k=k+1;
			count=0;
			fnzi[k]=i;
		}
	}
	fnzi[0]=0;
	sizeindi[procs-1]=A.nrows-sizeindi[procs-2]+2;			
}
	
    	MPI_Scatter(ncols, 1, MPI_INT, &ncolsm, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(nrows, 1, MPI_INT, &nrowsm, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(nozeros, 1, MPI_INT, &nozerosm, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(fnz, 1, MPI_INT, &fnzm, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(sizeindi, 1, MPI_INT, &sizeindim, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(fsearch, 1, MPI_INT, &fsearchm, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(lsearch, 1, MPI_INT, &lsearchm, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(sumnnz, 1, MPI_INT, &sumnnzm, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vector bm(nrowsm);
	matrix Am(nrowsm,ncolsm,nozerosm);
	
	MPI_Scatter(b.val, size/procs, MPI_DOUBLE, bm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Scatterv(A.val, nozeros, fnz, MPI_DOUBLE, Am.val, nozerosm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(A.indj, nozeros, fnz, MPI_INT, Am.indj, nozerosm, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Scatterv(A.indi, sizeindi, fnzi, MPI_INT, Am.indi, sizeindim, MPI_INT, 0, MPI_COMM_WORLD);		
	

	//Jacobi
	matrix Km(nrowsm,ncolsm,nrowsm);
	
	k=0;
	for(int i=0; i<(nrowsm+1); i++){
		for (int j=(Am.indi[i]-sumnnzm); j<(Am.indi[i+1]-sumnnzm); j++){
			if(Am.indj[j]==fsearchm){
				Km.val[k]=Am.val[j];
				Km.indj[k]=Am.indj[j];
				Km.indi[k]=Am.indj[j];
				k=k+1;
				fsearchm=fsearchm+1;
			
				break;
			}
		}
		
	}
	Km.indi[nrowsm]=(rank+1)*nrowsm;


	//bicgstab
	vector xm(bm.length),pm(bm.length),vm(bm.length),helpvm(bm.length),helpvm2(bm.length),rm(bm.length),rhatm(bm.length),ym(bm.length);
	vector sm(bm.length),zm(bm.length),tm(bm.length);
	vector x(b.length),p(b.length),v(b.length),helpv(b.length),helpv2(b.length),r(b.length),rhat(b.length),y(b.length),s(b.length);
	vector z(b.length),t(b.length);
	matrix Kinv(Km.nrows,Km.ncols,Km.nnz);
	double Nmax=1000, tol=pow(10.0,-10.0);
	double rho=1;
	double alpha=1;
	double omega=1;
	double rho2;
	double beta;
	double w1;
	double w2;
	double help;
	double help2;
	
	//arxikopoiisi x=0 v=0 p=0
	for (int i=0; i<bm.length; i++) { 
		xm.val[i]=0;
		vm.val[i]=0;
		pm.val[i]=0;
	} 

	//anastrofos tou Km dhladh tou K poy yparxei se ka8e procs
	inverse(Kinv,Km); 
	
	//residual
	MPI_Gather(xm.val, size/procs, MPI_DOUBLE, x.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	SpMV(helpvm,Am,x,sumnnzm);
	VVdic(rm,bm,helpvm);
	
	
	MPI_Gather(rm.val, size/procs, MPI_DOUBLE, r.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(r.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (int i=0; i<b.length; i++){
		rhat.val[i]=r.val[i];
		
	}	
	
	for (int iter=0; iter<100; iter++){
		rho2=rho;
		VVmul(rho,rhat,r);
		beta=( (rho/rho2)*(alpha/omega));

		CVmul(helpvm,omega,vm);
		VVdic(helpvm2,pm,helpvm);
		CVmul(helpvm,beta,helpvm2);
		MPI_Scatter(r.val, size/procs, MPI_DOUBLE, rm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		VVadd(pm,rm,helpvm);
		MPI_Gather(pm.val, size/procs, MPI_DOUBLE, p.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(p.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		SpMV(ym,Kinv,p,Kinv.indi[0]);
		
		MPI_Gather(ym.val, size/procs, MPI_DOUBLE, y.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(y.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		SpMV(vm,Am,y,sumnnzm);
		
		MPI_Gather(vm.val, size/procs, MPI_DOUBLE, v.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(v.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		VVmul(help,rhat,v);
		alpha=rho/help;

		MPI_Scatter(v.val, size/procs, MPI_DOUBLE, vm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		CVmul(helpvm,alpha,vm);
		

		VVdic(sm,rm,helpvm);
		MPI_Gather(sm.val, size/procs, MPI_DOUBLE, s.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(s.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		SpMV(zm,Kinv,s,Kinv.indi[0]);
		MPI_Gather(zm.val, size/procs, MPI_DOUBLE, z.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(z.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		SpMV(tm,Am,z,sumnnzm);

		MPI_Gather(tm.val, size/procs, MPI_DOUBLE, t.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(t.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		SpMV(helpvm,Kinv,t,Kinv.indi[0]);
		SpMV(helpvm2,Kinv,s,Kinv.indi[0]);
		MPI_Gather(helpvm.val, size/procs, MPI_DOUBLE, helpv.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(helpv.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(helpvm2.val, size/procs, MPI_DOUBLE, helpv2.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(helpv2.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		

		VVmul(w1,helpv,helpv2);
		VVmul(w2,helpv,helpv);
		omega=w1/w2;
		
		MPI_Scatter(y.val, size/procs, MPI_DOUBLE, ym.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(z.val, size/procs, MPI_DOUBLE, zm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(x.val, size/procs, MPI_DOUBLE, xm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		CVmul(helpvm,alpha,ym);
		CVmul(helpvm2,omega,zm);
		VVadd(helpvm,xm,helpvm);
		VVadd(xm,helpvm,helpvm2);
		MPI_Gather(xm.val, size/procs, MPI_DOUBLE, x.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(x.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Scatter(t.val, size/procs, MPI_DOUBLE, tm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(s.val, size/procs, MPI_DOUBLE, sm.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		CVmul(helpvm,omega,tm);
		VVdic(rm,sm,helpvm);

		MPI_Gather(rm.val, size/procs, MPI_DOUBLE, r.val, size/procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(r.val, b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		if (rank==0){
			for (int i=0; i<x.length; i++){
				cout<<x.val[i]<<endl;
		
			}	
			cout<<endl;
			norm(help,r);
			norm(help2,rhat);
		
			if (help<(tol*help2)){
				cout<<"Number of iterations is:"<<iter+1<<endl;
				cout<<"Time "<<time<<endl;
				break;
			}
		}
	}

	MPI_Finalize();
	return 0;
    	
}




#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

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

void SpMV(vector &y,const matrix &M, const vector &x){
	double sum=0;
	for(int i = 0 ; i<M.nrows;i++){
		sum=0;
		for (int j=M.indi[i];j<M.indi[i+1];j++){
			sum+=M.val[j]*x.val[M.indj[j]];
		}
		y.val[i]=sum;
	}
}
void SpMadd(matrix &C,matrix &A,matrix &B){
	int w[A.nrows],tnnz=0;
	double ww[A.nrows];
	int ldg=-2;
	for(int i=0;i<A.nrows;i++){
		w[i]=-1;
		ww[i]=0;
	}
	for(int i=0;i<A.nrows;i++){
		for(int j=A.indi[i];j<A.indi[i+1];j++){
			w[A.indj[j]]=ldg;
			ldg=A.indj[j];
			tnnz++;
		}
		for(int j=B.indi[i];j<B.indi[i+1];j++){
			if(w[B.indj[j]]!=-1){
				w[B.indj[j]]=ldg;
				ldg=B.indj[j];
				tnnz++;
			}
		}
	}
	C.nnz=tnnz;
	C.val = new double[C.nnz];
	C.indj = new int[C.nnz];
	for(int i=0;i<A.nrows;i++){
		for(int j=A.indi[i];j<A.indi[i+1];j++){
			w[A.indj[j]]=ldg;
			ldg=A.indj[j];
			tnnz++;
		}
		for(int j=B.indi[i];j<B.indi[i+1];j++){
			if(w[B.indj[j]]!=-1){
				w[B.indj[j]]=ldg;
				ldg=B.indj[j];
				tnnz++;
			}
		}
	}
}

void VVadd(vector &y,const vector &a, const vector &b){
	for(int i = 0 ; i<y.length; i++){
		y.val[i]=a.val[i]+b.val[i];
	}
}

void VVdic(vector &y,const vector &a, const vector &b){
	for(int i = 0 ; i<y.length; i++){
		y.val[i]=a.val[i]-b.val[i];
	}
}

void VVmul(double &y,const vector &a, const vector &b){
        double sum=0;
	for(int i = 0 ; i<a.length; i++){
		sum=sum+a.val[i]*b.val[i];
	}
	y=sum;
}

void CVmul(vector &y,double &a, const vector &b){
	for(int i = 0 ; i<y.length; i++){
		y.val[i]=a*b.val[i];
	}
}

void inverse(matrix &K,const matrix &A){
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

void Bicgstab(vector &x,const matrix &A,const vector &b,const matrix &K,double &tol,int &Nmax){
	vector r(b.length);
	vector rhat(b.length);
	vector v(b.length);
	vector p(b.length);
	vector y(b.length);
	vector z(b.length);
	vector t(b.length);
	vector s(b.length);
	vector helpv(b.length);
	vector helpv2(b.length);
	matrix Kinv(K.nrows,K.ncols,K.nnz);
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
	for (int i=0; i<b.length; i++) { 
		x.val[i]=0;
		v.val[i]=0;
		p.val[i]=0;
		//cout<<"sto diaolo";
	} 
	//antistrofi toy K. tha xreiastei stin poreia
	inverse(Kinv,K);
	
	//residual 
	SpMV(helpv,A,x);
	VVdic(r,b,helpv);
	for (int i=0; i<b.length; i++){
		rhat.val[i]=r.val[i];
		
	}
	//epilysi
	for (int iter=0; iter<Nmax; iter++){
		rho2=rho;
		VVmul(rho,rhat,r);
		beta=( (rho/rho2)*(alpha/omega));
		
		CVmul(helpv,omega,v);
		VVdic(helpv2,p,helpv);
		CVmul(helpv,beta,helpv2);
		VVadd(p,r,helpv);

		SpMV(y,Kinv,p);
		SpMV(v,A,y);
		
		VVmul(help,rhat,v);
		alpha=rho/help;
		
		CVmul(helpv,alpha,v);
		VVdic(s,r,helpv);
		
		SpMV(z,Kinv,s);
		SpMV(t,A,z);
		
		SpMV(helpv,Kinv,t);
		SpMV(helpv2,Kinv,s);
		VVmul(w1,helpv,helpv2);
		VVmul(w2,helpv,helpv);
		omega=w1/w2;
		
		
		CVmul(helpv,alpha,y);
		CVmul(helpv2,omega,z);
		VVadd(helpv,x,helpv);
		VVadd(x,helpv,helpv2);
		
		
		CVmul(helpv,omega,t);
		VVdic(r,s,helpv);
		
		for (int i=0; i<x.length; i++){
			cout<<x.val[i]<<endl;
		
		}	
		cout<<endl;


		norm(help,r);
		
		norm(help2,rhat);
		
		if (help<(tol*help2)){
			cout<<"Number of iterations is:"<<iter+1<<endl;
			break;
		}
	}
	
	

}

int main(void){
	
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
    	matrix C(size,size,nnz); // dhlwnoume ton pinaka 
    	getline(file, str);   // trith grammh gia to val[]  
    	std::stringstream stream3(str);   
    	double v;  
    	int i=0;      
    	while(stream3>>v){                
         	if(!stream) break;     
         	C.val[i]=v;
         	i++;                                          
    	}  
    	getline(file, str);   // 4h grammh gia to indj[]  
    	std::stringstream stream4(str);        
    	i=0;      
    	while(stream4>>n){                
         	if(!stream) break;     
         	C.indj[i]=n;
         	i++;                                          
    	} 
    	getline(file, str);   // 5h grammh gia to indi[]  
    	std::stringstream stream5(str);        
    	i=0;      
    	while(stream5>>n){                
         	if(!stream) break;     
         	C.indi[i]=n;
         	i++;                                          
    	}         
	std::ifstream file2("b.txt"); // dhmiourgoume to vector 
    	std::string str2; 
    	getline(file2, str2);
    	std::stringstream streamB1(str2);
    	streamB1>>n;
    	size=n;
    	vector z(size);
    	getline(file2, str2);   // fortonoume tis times sto vector mas   
    	std::stringstream streamB2(str2);        
    	i=0;      
    	while(streamB2>>v){                
         	if(!stream) break;     
         	z.val[i]=v;
         	i++;                                          
    	} 
    
    	
	//prakseis
    	matrix M(size, size,size);
	vector x(size);
	double tol=pow(10.0,-10.0);
	int max=80;
	

	Jacobi(M,C);
	Bicgstab(x,C,z,M,tol,max);

	
    	
}





















#define N 7
#define NRATE 20
#include<math.h>
#include<stdio.h>
#include<iostream>
#include<string.h>
using namespace std;
extern double aion[];
extern double zion[];
enum elements {he4,c12,o16,ne20,mg24,si28,ni56};
enum rates {ircag,iroga,ir3a,irg3a,ir1212,ir1216,ir1616,iroag,irnega,irneag,irmgga,irmgag,irsiga,ircaag,irtiga,irni2si,irsi2ni,irsi2nida,irni2sida,irsi2nidsi};   
class Linear_system{
	private :
	double a[N][N];
	int nzr[N];
	int	nzc[N];
	int ri[N][N], ci[N][N];
	bool nz[N][N];
	int piv[N];
	double b[N];
	
	public :
	Linear_system(){
		for(int row=0;row<N;row++)
			for(int col=0;col<N;col++){
				a[row][col]=0;
				ri[row][col]=0;
				ci[row][col]=0;
				nz[row][col]=0;
			}
		for(int ind=0;ind<N;ind++){
			nzr[ind]=0;
			nzc[ind]=0;
			piv[ind]=0;
			b[ind]=0;
		}
//		for(int row=0;row<N;row++)
//			for(int col=0;col<N;col++)
//				if(A[row][col]!=0) *GetElementAddress(row,col)=A[row][col];
		
	}
	void initial(double A[N][N]){
	
		for(int row=0;row<N;row++)
			for(int col=0;col<N;col++)
				if(A[row][col]!=0) *GetElementAddress(row,col)=A[row][col];
//				a[row][col]=ap[row][col];
//		for(int ind=0;ind<N;ind++)
//			b[ind]=B[ind];
	}
		
	double* GetElementAddress(int row, int col) {
	if (!(nz[row][col])) {
	nz[row][col] = 1;
	ri[col][nzr[col]] = row;
	ci[row][nzc[row]] = col;
	nzr[col]++;
	nzc[row]++;
	a[row][col] = 0.;
	}
	return a[row] + col;
	}
	                                          
	void CleanMatrix() {
	int i, k;
	for (i = 0; i < N; i++) {
		for (k = 0; k < N; k++) {
			nz[i][k] = 0;
		}
	}
	memset(nzr, 0, sizeof(int)*N);
	memset(nzc, 0, sizeof(int)*N);
	}
	

	int FindPiv(int i, char *todo, double minpiv) {
		int k, pvr, nc;
		double mul, mulpiv, fari;
		nc = N+1;
		pvr = ri[i][0];
		mul = 0.;
		for (k = 0; k < nzr[i]; k++) {
			if (todo[ri[i][k]]) {
			fari = fabs(a[ri[i][k]][i]);
			if (fari > mul) {
			mul = fari;
			}
		}
	}
	mulpiv = mul*minpiv;
	for (k = 0; k < nzr[i]; k++) {
		if (todo[ri[i][k]]) {
		fari = fabs(a[ri[i][k]][i]);
			if (fari >= mulpiv && nzc[ri[i][k]] < nc) {
			nc = nzc[ri[i][k]];
			pvr = ri[i][k];
			}
		}
	}
	return pvr;
	}
	
	void naivfct() {
		char todo[N];
		int i, j, k, pvr;
		double den, mul;
		for (k = 0; k < N; k++) {
			todo[k] = 1;
		}
		for (i = 0; i < N; i++) {
			piv[i] = pvr = FindPiv(i, todo, 1.E-8);
			todo[pvr] = 0;
			den = 1.0/a[pvr][i];
			a[pvr][i] = den;
		for (k = 0; k < nzr[i]; k++) {
			if (!todo[ri[i][k]]) { continue; }
			mul = -a[ri[i][k]][i]*den;
			a[ri[i][k]][i] = mul;
		for (j = 0; j < nzc[pvr]; j++) {
			if (ci[pvr][j] <= i) { continue; }
			(*GetElementAddress(ri[i][k],ci[pvr][j]))+=mul*a[pvr][ci[pvr][j]];
	//		a[ri[i][k]][ci[pvr][j]]+=mul*a[pvr][ci[pvr][j]];
		}
		}
		}
//		invpiv={3,4,0,1,2,5,6};
	}
	void backward(double B[N],double x[N]){
		
		for(int ind=0;ind<N;ind++)
			b[ind]=B[ind];
	
		for(int k=0;k<N;k++){
			for(int i=k+1;i<N;i++){
				int pvr=piv[i];
				if(nz[pvr][k])
					b[pvr]+=a[pvr][k]*b[piv[k]];
			}
		}
	
		b[piv[N-1]]*=a[piv[N-1]][N-1];
		for(int i=N-2;i>=0;i--){
			double total=0;
			int pvr=piv[i];
			for(int j=i+1;j<N;j++)
				total+=a[pvr][j]*b[piv[j]];
			b[pvr]=(b[pvr]-total)*a[pvr][i];
		}
		for(int i=0;i<N;i++)
			x[i]=b[piv[i]];
	}	
	void print(double ap[N][N]){
		for(int row=0;row<N;row++){
			for(int col=0;col<N;col++)
				cout<<a[row][col]<<"   ";
			cout<<"\n";
		}
		for(int ind=0;ind<N;ind++)
			cout<<piv[ind]<<"     ";
		cout<<"\n";
		for(int ind=0;ind<N;ind++)
			cout<<b[ind]<<"\n";
		cout<<endl;
		for(int i=0;i<N;i++){
			double total=0;
			for(int j=0;j<N;j++)
				total+=b[piv[j]]*ap[i][j];
			cout<<total<<endl;
		}
	}
	void add(double y[N]){
		for(int i=0;i<N;i++)
			y[i]+=b[piv[i]];
	}
};
void screen(double,double,double [],double[],double []);
void get_rate(double,double ,double [],double []);
void jacob(double [],double [],double [N][N]);
void ydot(double [],double [],double []);
void BD(const double ,const double ,double [],const double);

#include "nuclear.h"
#include "../paul.h"
extern double RHO_CONST;
extern double P_CONST;
#define TEMP_CONST (P_CONST/RHO_CONST)
const double C=299792458;

double realmass[N]={4.002602,12.0,15.9949,19.9924401754,23.985041700,27.9769265325,55.942132};
extern "C"{
	void burning(struct domain *,double);
	double  get_dV(double *,double *);
	}
void cell_burning(const double den,const double temp,double *x,const double t_end,double &average_rate){
        double y[N];
	double oldmass=0;
        for(int i=0;i<N;i++){
                y[i]=x[i]/aion[i];
		oldmass+=y[i]*realmass[i];
	}
        double enr=0;
	double newmass=0;
        BD(den,temp,y,t_end);
        for(int i=0;i<N;i++){
                x[i]=y[i]*aion[i];
		newmass+=y[i]*realmass[i];
//		cout<<x[i]<<"   ";
	}
//	cout<<endl;

	enr=(oldmass-newmass)*C*C*1e4;
	average_rate=enr/t_end;
		
        }
void burning(struct domain *theDomain,double dt){
	int Nr = theDomain->Nr;
	int Nz = theDomain->Nz;
	int * Np = theDomain->Np;
	int Npl = theDomain->Npl;
	struct cell ** theCells = theDomain->theCells;
	double * r_jph = theDomain->r_jph;
   	double * z_kph = theDomain->z_kph;
	
	int Ng = theDomain->Ng;
	int * dim_rank = theDomain->dim_rank;
        int * dim_size = theDomain->dim_size;
        int Z_Periodic = theDomain->theParList.Z_Periodic;


   	int i,j,k,p;

        int j_min = 0;
   	int j_max = Nr;
   	int k_min = 0;
   	int k_max = Nz;
	if( dim_rank[0] != 0 ) j_min = Ng;
   	j_max = Nr-Ng;
	if( dim_rank[1] != 0 || Z_Periodic ) k_min = Ng;
        if( dim_rank[1] != dim_size[1]-1 || Z_Periodic ) k_max = Nz-Ng;


   	for( j=j_min ; j<j_max ; ++j ){
      		for( k=0 ; k<Nz ; ++k ){
         		int jk = j+Nr*k;
         		for( i=0 ; i<Np[jk] ; ++i ){
            		struct cell * c = &(theCells[jk][i]);
            		double phip = c->piph;
            		double phim = c->piph - c->dphi;
            		double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            		double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            		double dV = get_dV(xp,xm);
			double xmass[N];
			for(int i=0;i<N;i++)
				xmass[i]=c->cons[NUM_C+i]/c->cons[0];			
//			for(int i=0;i<N;i++)
//				cout<<c->prim[NUM_C+i]<<"   ";
//			cout<<endl;
//			
//			if(xp[0]<7e-2){
//				for(int i=0;i<N;i++)
//					cout<<xmass[i]<<"   ";
//				cout<<endl;
//			}
         		double average_rate=0;
         		double temp=c->temp;
         		double rho=c->prim[RHO];
         		cell_burning(rho*RHO_CONST,temp,xmass,dt,average_rate);
			for(int i=0;i<N;i++)
				c->cons[NUM_C+i]=xmass[i]*c->cons[0];
//			cout<<endl;
//			cout<<dt<<"   "<<temp<<"   "<<c->cons[TAU]<<"   "<<average_rate*dt/TEMP_CONST*rho*dV<<endl;
         		c->cons[TAU]+=average_rate*dt/TEMP_CONST*rho*dV;

//			if(xp[0]<7e-2){
//                                for(int i=0;i<N;i++)
//                                        cout<<xmass[i]<<"   ";
//                                cout<<endl;
//                        }

    }
   }
  }
}
/*
double mindt_burn(double *prim,double temp,double enr){
        double energy_rate;
        double rho=prim[RHO];
        double comp[2]; 
        comp[he4]=prim[NUM_C+he4];
        comp[c12]=prim[NUM_C+c12];
        cell_burning(rho*RHO_COF,temp*TEMP_COF,comp,0,&energy_rate);
//      printf("%f   %f   %f\n",rho,temp,enr);
        return fabs(0.1*enr/(energy_rate/TEMP_COF));
//      return 100;
}
*/

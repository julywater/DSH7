//#include <iostream>

#include "../paul.h"
#include "../helm.h"


double RHO_CONST=1e6;
double P_CONST=1e23;
#define TEMP_CONST (P_CONST/RHO_CONST)
#define VEL_CONST (sqrt(P_CONST/RHO_CONST))
extern double * table;

//using namespace std;

// extern "C" {
     
//      void teos_(double *,double *,double *,double *,double *,double *,double *,double *);

//        }
double min(double a,double b){
	if(a<=b)
		return a;
	else 
		return b;
}
double max(double a,double b){
	if(a>=b)
		return a;
	else 
		return b;
}
double get_temp(double *prim){
	return prim[PPP]/prim[RHO]*1e9;
}
void invert_helmeos_ed(double ein,double dens,double *temp,double *pres,double *cs){
	ein=ein/dens*TEMP_CONST; //energy density per mass
        dens=dens*RHO_CONST;
	double eostol=1e-8;
	double slmin=1e-14;
	double ein_temp;
	double det,dpt;

	teos_(temp,&dens,&ein_temp,pres,cs,&det,&dpt,table);
	
	double f=ein_temp/ein-1.0;
        double df=det/ein;
	double eos_slope=f/df;
        double tmpnew = min(max(0.5*(*temp),*temp-eos_slope),2.0*(*temp));
	double tmpdif=fabs((tmpnew - *temp)/(*temp));
        *temp=min(1.0e14,max(tmpnew,1.0e-11));
	
 	int i=0;
	while(tmpdif>eostol&&fabs(eos_slope)>slmin&&i<40){
        
		teos_(temp,&dens,&ein_temp,pres,cs,&det,&dpt,table);
                f=ein_temp/ein-1.0;
       		df=det/ein;
		eos_slope=f/df;
		tmpnew = min(max(0.5*(*temp),*temp-eos_slope),2.0*(*temp));
		tmpdif=fabs((tmpnew - *temp)/(*temp));
		*temp=min(1.0e14,max(tmpnew,1.0e-11));
		i++;
	}

//	teos_(temp,&dens,&ein_temp,pres,cs,&det,&dpt,table);
        *cs=*cs/VEL_CONST;
        *pres=*pres/P_CONST;

	}


void invert_helmeos_pd(double pin,double dens,double *temp,double *eout,double *cs){
        pin=pin*P_CONST;
        dens=dens*RHO_CONST;
	double eostol=1e-8;
	double slmin=1e-14;
	double p_temp=0;
	double det,dpt;
	
//	printf("%f   %f   %f\n",dens,pin,*temp);
	teos_(temp,&dens,eout,&p_temp,cs,&det,&dpt,table);
	

	double f=p_temp/pin-1.0;
	double df=dpt/pin;
	double eos_slope=f/df;
	double tmpnew = min(max(0.5*(*temp),*temp-eos_slope),2.0*(*temp));
	double tmpdif=fabs((tmpnew - *temp)/(*temp));
        *temp=min(1.0e14,max(tmpnew,1.0e-11));

	int i=0;
	while(tmpdif>eostol&&fabs(eos_slope)>slmin&&i<40){
        
		teos_(temp,&dens,eout,&p_temp,cs,&det,&dpt,table);
                f=p_temp/pin-1.0;
       		df=dpt/pin;
		eos_slope=f/df;
		tmpnew = min(max(0.5*(*temp),*temp-eos_slope),2.0*(*temp));
		tmpdif=fabs((tmpnew - *temp)/(*temp));
		*temp=min(1.0e14,max(tmpnew,1.0e-11));
		i++;
	}

//	teos_(temp,&dens,eout,&p_temp,cs,&det,&dpt,table);
	
        *eout=*eout/TEMP_CONST*(dens/RHO_CONST);

	*cs=*cs/VEL_CONST;
	
	}

	

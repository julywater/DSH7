#include<iostream>
#include<math.h>
#include<time.h>
#include<stdio.h>
#define N 7
#define NRATE 20
using namespace std;
enum elements {he4,c12,o16,ne20,mg24,si28,ni56};
enum rates {ircag,iroga,ir3a,irg3a,ir1212,ir1216,ir1616,iroag,irnega,irneag,irmgga,irmgag,irsiga,ircaag,irtiga,irni2si,irsi2ni,irsi2nida,irni2sida,irsi2nidsi};   
double zion[N]={2,6,8,10,12,14,28};
double aion[N]={4,12,16,20,24,28,56};
void screen5(const double,const double,const double,const double,const double,double [],double [],double [],double [],double [],double[],double[]);

void get_rate(double temp,double den,double y[N],double rate[NRATE]){
	int i;
	double tt9,t9r,t9,t912,t913,t923,t943,t953,t932,t92,t93,t972,t9r32,t9i,t9i13,t9i23,t9i32,t9i12,t9ri,term,term1,term2,term3,rev,r2abe,t9a,t9a13,t9a56,t9a23,gt9h,rbeac,oneth,fivsix;
	oneth = 1.0/3.0;
	fivsix=5.0/6.0;
	for(int ind=0;ind<NRATE;ind++)
		rate[ind]=0;
	
//some temperature factors
//limit t9 to 10, except for the reverse ratios
	tt9 = temp * 1.0e-9;
	t9r = tt9;
	t9 = min(tt9,10.0);
	t912 = sqrt(t9);
	t913 = pow(t9,oneth);
	t923 = t913*t913;
	t943 = t9*t913;
	t953 = t9*t923;
	t932 = t9*t912;
	t92 = t9*t9;
	t93 = t9*t92;
	t972 = t912*t93;
	t9r32 = t9r*sqrt(t9r);
	t9i = 1.0/t9;
	t9i13 = 1.0/t913;
	t9i23 = 1.0/t923;
	t9i32 = 1.0/t932;
	t9i12 = 1.0/t912;
	t9ri = 1.0/t9r;
   
   
//triple alpha to c12
	r2abe=(7.40e5*t9i32)*exp(-1.0663*t9i)+4.164e9*t9i23*exp(-13.49*t9i13 - t92/0.009604)*(1.0+0.031*t913+8.009*t923+1.732*t9+49.883*t943+27.426*t953);
	rbeac=(130.*t9i32)*exp(-3.3364*t9i)+2.510e7*t9i23*exp(-23.57*t9i13-t92/0.055225)*(1.0+0.018*t913+5.249*t923+0.650*t9 +19.176*t943+6.034*t953);
	if (t9>0.08) 
		term = 2.90e-16 * (r2abe*rbeac)+0.1 * 1.35e-07*t9i32*exp(-24.811*t9i);
	else
		term = 2.90e-16*(r2abe*rbeac)*(0.01+0.2*(1.0+4.0*exp(-pow((0.025*t9i),3.263)))/(1.0 + 4.0*exp(-pow((t9/0.025),9.227))))+0.1 * 1.35e-07 * t9i32 * exp(-24.811*t9i);

	rate[ir3a]=term*(den*den)/6.0;
	rev=2.00e20*exp(-84.424*t9ri);
	rate[irg3a]=rev*(t9r*t9r*t9r)*term;
	
//c12 + c12 reaction from caughlan and fowler 1988
	t9a = t9/(1.0+0.0396*t9);
	t9a13 = pow(t9a,oneth);
	t9a56 = pow(t9a,fivsix);
	term = 4.27e+26 * t9a56/t932*exp(-84.165/t9a13-2.12e-03*t9*t9*t9);
	rate[ir1212]= 0.5*den*term;
	
//c12 + o16 reaction
	if (t9>0.5) {
		t9a = t9/(1.+0.055*t9);
		t9a13 = pow(t9a,oneth);
		t9a23 = t9a13*t9a13;
		t9a56 = pow(t9a,fivsix);
		term = 1.72e31 * t9a56 * t9i32 * exp(-106.594/t9a13) /(exp(-0.18*t9a*t9a) + 1.06e-03*exp(2.562*t9a23));
		rate[ir1216]=den*term;
	}
	else
		rate[ir1216]= 0.0;
		
//16o+16o rate
	term=7.10e36 * t9i23*exp(-135.93 * t9i13 - 0.629*t923-0.445*t943 + 0.0103*t9*t9);
	rate[ir1616]=0.5*den*term;
	
//12c(ag)16o and inverse
	term = 1.04e8/(t92*pow((1.0+ 0.0489*t9i23),2))*exp(-32.120*t9i13-t92/12.222016)+1.76e8/(t92*pow((1.0+0.2654*t9i23),2))*exp(-32.120*t9i13)+1.25e3*t9i32*exp(-27.499*t9i)+ 1.43e-2*t92*t93*exp(-15.541*t9i);
	term = 1.7* term;
	rate[ircag]=term * den;
	rev = 5.13e10 * t9r32 * exp(-83.111*t9ri);
	rate[iroga]= rev * term;

	
//16o(ag)20ne + inverse
	term1 = 9.37e9*t9i23*exp(-39.757*t9i13 - t92/2.515396);
	term2 = 62.1 * t9i32 * exp(-10.297*t9i) +538.0 * t9i32 * exp(-12.226*t9i) +13* t92 * exp(-20.093*t9i);
	term=term1 + term2;
	rate[iroag]= den * term;
	rev=5.65e10*t9r32*exp(-54.937*t9ri);
	rate[irnega]= rev * term;
	
//20ne(ag)24mg + inverse
	term1 = 4.11e11*t9i23*exp(-46.766*t9i13 - t92/4.923961)*(1.0+0.009*t913+0.882*t923+0.055*t9+0.749*t943 + 0.119*t953);
	term2 = 5.27e03*t9i32*exp(-15.869*t9i)+6.51e03*t912*exp(-16.223*t9i);
	term3 = 0.1*(42.1*t9i32*exp(-9.115*t9i)+32.0*t9i23*exp(-9.383*t9i));
	term = (term1+term2+term3)/(1.0+5.0*exp(-18.960*t9i));
	rate[irneag]=den*term;
	rev=6.01e10*t9r32*exp(-108.059*t9ri);
	rate[irmgga]= term * rev;

//24mg(ag)28si + inverse
	term1=(1.0+5.0*exp(-15.882*t9i));
	term=(4.78e1*t9i32*exp(-13.506*t9i)+2.38e+03*t9i32*exp(-15.218*t9i)+2.47e+02*t932*exp(-15.147*t9i)+0.1*(1.72e-9*t9i32 * exp(-5.028*t9i)+1.25e-03 * t9i32 * exp(-7.929*t9i)+2.43e01 * t9i * exp(-11.523*t9i)))/term1;
	rate[irmgag]= den * term;
	rev = 6.27e10 * t9r32 * exp(-115.862*t9ri);
	rate[irsiga]=rev * term;
	
//c..40ca(ag)44ti + inverse
	term=4.66e+24*t9i23*exp(-76.435 * t9i13 *(1.0+ 1.650e-02*t9 + 5.973e-03*t92 -3.889e-04*t93));
	rev=6.843e+10*t9r32*exp(-59.510*t9ri);
	rate[ircaag]= den*term;
	rate[irtiga]= rev*term;
	
//
//	

	if (t9>2.5){
		rate[irsi2ni]=(pow(t9i32,3))*exp(239.42*t9i-74.741)*pow(den,3)*rate[ircaag];
		rate[irsi2nida]=3.0*rate[irsi2ni];
		rate[irsi2nidsi]=rate[irsi2ni];
		rate[irni2si]= pow(t932,3)* exp(-274.12*t9i+74.914) * rate[irtiga]/pow(den,3);
		rate[irni2sida]=-3.0*rate[irni2si];
		
	}


//..the user may want to insert screening factors for the relevant rates here
}
void screen(double temp,double den,double y[],double rate0[],double rate[]){
	for(int i=0;i<NRATE;i++)
		rate[i]=rate0[i];
	double abar=0;
	double zbar=0;
	double z2bar=0;
	for(int i=0;i<7;i++){
       abar+=y[i];
       zbar += zion[i]*y[i];
       z2bar+= zion[i]*zion[i]*y[i];
    }
    abar=1.0/abar;
    zbar= zbar*abar;
    z2bar=z2bar*abar;
	double sca[9],scadt[9],scadd[9];
	double z1[9]={zion[he4],zion[he4],zion[c12],zion[c12],zion[c12],zion[o16],zion[o16],zion[ne20],zion[mg24]};
	double a1[9]={aion[he4],aion[he4],aion[c12],aion[c12],aion[c12],aion[o16],aion[o16],aion[ne20],aion[mg24]};
	double z2[9]={zion[he4],4,zion[he4],zion[c12],zion[o16],zion[o16],zion[he4],zion[he4],zion[he4]};
	double a2[9]={aion[he4],4,aion[he4],aion[c12],aion[o16],aion[o16],aion[he4],aion[he4],aion[he4]};
	screen5(temp,den,zbar,abar,z2bar,z1,a1,z2,a2,sca,scadt,scadd);
	rate[ir3a]*=sca[0]*sca[1];
	rate[ircag]*=sca[2];
	rate[ir1212]*=sca[3];
	rate[ir1216]*=sca[4];
	rate[ir1616]*=sca[5];
//	rate[ir1216] possible bug in Timmes
	rate[iroag]*= sca[6];
	rate[irneag]*=sca[7];
	rate[irmgag]*=sca[8];
/*
	if(y[c12]+y[o16]>0.004) {
		rate[irsi2ni]=0;
       		rate[irsi2nida]=0;
        	rate[irsi2nidsi]=0;
      		rate[irni2si]=0;
        	rate[irni2sida]=0;
	}
*/		
/*
	double sc1a,sc2a,scadt,scadd;
	screen5(temp,den,zbar,abar,z2bar,zion[he4],aion[he4],zion[he4],aion[he4],0,0,sc1a,scadt,scadd);
    screen5(temp,den,zbar,abar,z2bar,zion[he4],aion[he4],4.0,8.0,0,0,sc2a,scadt,scadd);
	rate[ir3a]*=sc1a*sc2a;
	screen5(temp,den,zbar,abar,z2bar, zion[c12],aion[c12],zion[he4],aion[he4],0,0,sc1a,scadt,scadd);
    rate[ircag]*=sc1a;
	screen5(temp,den,zbar,abar,z2bar,zion[c12],aion[c12],zion[c12],aion[c12], 0,0,sc1a,scadt,scadd);
	rate[ir1212]*=sc1a;
	screen5(temp,den,zbar,abar,z2bar,zion[c12],aion[c12],zion[o16],aion[o16],0,0,sc1a,scadt,scadd);
	rate[ir1216]*=sc1a;
	screen5(temp,den,zbar,abar,z2bar,zion[c12],aion[c12],zion[o16],aion[o16],0,0,sc1a,scadt,scadd);
	rate[ir1216] *=sc1a;
	screen5(temp,den,zbar,abar,z2bar,zion[o16],aion[o16],zion[he4],aion[he4],0,0,sc1a,scadt,scadd);
	ratdum(iroag)    = ratraw(iroag) * sc1a
	screen5(temp,den,zbar,abar,z2bar,zion[ne20],aion[ne20],zion[he4],aion[he4],0,0,sc1a,scadt,scadd);
    rate[irneag]*=sc1a;
	screen5(temp,den,zbar,abar,z2bar,zion[mg24],aion[mg24],zion[he4],aion[he4],0,0,sc1a,scadt,scadd);
    rate[irmgag]*=sc1a;
	
*/
	}
void jacob(double rate[NRATE],double y[N],double dfdy[N][N]){

//	
	double  rate_irsi2ni=0;
        double  rate_irsi2nida=0;
        double  rate_irsi2nidsi=0;
        double  rate_irni2si=0;
        double  rate_irni2sida=0;

	if ((y[c12]+y[o16])<0.004){
		rate_irsi2ni=rate[irsi2ni]*pow(y[he4],3)*y[si28];
		rate_irsi2nida=rate[irsi2nida]*pow(y[he4],2)*y[si28];
		rate_irsi2nidsi=rate[irsi2nidsi]*pow(y[he4],3);
		rate_irni2si=min(rate[irni2si]*pow(y[he4],3),1e20);
		if(rate_irni2si==1e20)
			rate_irni2sida=0;
		else 
			rate_irni2sida=rate[irni2sida]*pow(y[he4],2);
	}

	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			dfdy[i][j]=0;
	dfdy[he4][he4] = -6.0* y[he4] * y[he4] * rate[ir3a]-y[c12]*rate[ircag]-y[o16]* rate[iroag]-y[ne20]*rate[irneag]- y[mg24]*rate[irmgag]-7.0* rate_irsi2ni- 7.0*rate_irsi2nida*y[he4]+7*rate_irni2sida*y[ni56];
	dfdy[he4][c12] = 3.0*rate[irg3a]-y[he4]*rate[ircag] + 2.0 * y[c12]*rate[ir1212]+0.5* y[o16]*rate[ir1216];
	dfdy[he4][o16] = rate[iroga]+0.5*y[c12]*rate[ir1216]+2.0*y[o16]*rate[ir1616]-y[he4]*rate[iroag];
	dfdy[he4][ne20] = rate[irnega]-y[he4]*rate[irneag];
	dfdy[he4][mg24] = rate[irmgga]-y[he4]*rate[irmgag];
	dfdy[he4][si28] = rate[irsiga] -7.0* rate_irsi2nidsi*y[he4];
	dfdy[he4][ni56] =7.0*rate_irni2si;

	dfdy[c12][he4]=3.0*y[he4]*y[he4]*rate[ir3a]-y[c12]*rate[ircag];
	dfdy[c12][c12]=-rate[irg3a]-y[he4]*rate[ircag]-4.0*y[c12]*rate[ir1212]-y[o16]* rate[ir1216];
	dfdy[c12][o16]=rate[iroga]-y[c12]*rate[ir1216];
	
	
	dfdy[o16][he4]=y[c12]*rate[ircag]-y[o16]*rate[iroag];
	dfdy[o16][c12]=y[he4]*rate[ircag]-y[o16]*rate[ir1216];
	dfdy[o16][o16]= -rate[iroga]-y[c12]*rate[ir1216]-4.0* y[o16]*rate[ir1616]-y[he4]*rate[iroag];
	dfdy[o16][ne20]=rate[irnega];
	
	
	dfdy[ne20][he4] = y[o16]*rate[iroag]-y[ne20]*rate[irneag];
	dfdy[ne20][c12] = 2.0*y[c12]*rate[ir1212];
	dfdy[ne20][o16] =y[he4]*rate[iroag];
	dfdy[ne20][ne20] =-rate[irnega]-y[he4]*rate[irneag];
	dfdy[ne20][mg24] = -rate[irmgga];
			
	
	dfdy[mg24][he4]= y[ne20]*rate[irneag]-y[mg24]*rate[irmgag];
	dfdy[mg24][c12]= 0.5* y[o16]*rate[ir1216];
	dfdy[mg24][o16] = 0.5* y[c12]*rate[ir1216];
	dfdy[mg24][ne20]= y[he4]*rate[irneag];
	dfdy[mg24][mg24]= -rate[irmgga]-y[he4]*rate[irmgag];
	dfdy[mg24][si28]= rate[irsiga];
	
	
	dfdy[si28][he4]= y[mg24]*rate[irmgag]-rate_irsi2ni-rate_irsi2nida*y[he4]+rate_irni2sida*y[ni56];
	dfdy[si28][c12]=0.5*y[o16]*rate[ir1216];
	dfdy[si28][o16]=2.0*y[o16]*rate[ir1616]+0.5*y[c12]*rate[ir1216];
	dfdy[si28][mg24]=y[he4]*rate[irmgag];
	dfdy[si28][si28]=-rate[irsiga]-rate_irsi2nidsi*y[he4];
	dfdy[si28][ni56]=rate_irni2si;
	
	dfdy[ni56][he4]=rate_irsi2ni+rate_irsi2nida*y[he4]-rate_irni2sida*y[ni56];
	dfdy[ni56][si28]=rate_irsi2nidsi*y[he4];
	dfdy[ni56][ni56]=-rate_irni2si;
}
void ydot(double rate[NRATE],double y[N],double dydt[N]){

	double  rate_irsi2ni=0;
        double  rate_irsi2nida=0;
        double  rate_irsi2nidsi=0;
        double  rate_irni2si=0;
        double  rate_irni2sida=0;

        if ((y[c12]+y[o16])<0.004){
                rate_irsi2ni=rate[irsi2ni]*pow(y[he4],3)*y[si28];
                rate_irsi2nida=rate[irsi2nida]*pow(y[he4],2)*y[si28];
                rate_irsi2nidsi=rate[irsi2nidsi]*pow(y[he4],3);
                rate_irni2si=min(rate[irni2si]*pow(y[he4],3),1e20);
                if(rate_irni2si==1e20)
                        rate_irni2sida=0;
                else
                        rate_irni2sida=rate[irni2sida]*pow(y[he4],2);
        }


	dydt[he4]= 3.0 * y[c12]* rate[irg3a]\
- 3.0 * y[he4]*y[he4]*y[he4]* rate[ir3a]+ y[o16]* rate[iroga]- y[c12]*y[he4]* rate[ircag]\
+ y[c12]*y[c12]*rate[ir1212]+0.5* y[c12]*y[o16]*rate[ir1216]+y[o16]*y[o16]*rate[ir1616]\
- y[o16]* y[he4]*rate[iroag]+y[ne20]*rate[irnega]+y[mg24]* rate[irmgga]-y[ne20]* y[he4]*rate[irneag]\
+ y[si28]*rate[irsiga]-y[mg24]*y[he4]*rate[irmgag]-7.0*rate_irsi2ni*y[he4]+7.0*rate_irni2si*y[ni56];

	
	dydt[c12]= y[he4]*y[he4]*y[he4]*rate[ir3a] - y[c12]* rate[irg3a] + y[o16] * rate[iroga] - y[c12]*y[he4]*rate[ircag]-2.0* y[c12]*y[c12]*rate[ir1212]\
- y[c12]*y[o16]*rate[ir1216];

	dydt[o16]= - y[o16]*rate[iroga]+y[c12]*y[he4]*rate[ircag]-y[c12]*y[o16]*rate[ir1216]-2.0*y[o16]*y[o16]*rate[ir1616]-y[o16]*y[he4]*rate[iroag]+y[ne20]*rate[irnega];

	dydt[ne20]= y[c12]*y[c12]*rate[ir1212]+y[o16]*y[he4]*rate[iroag]-y[ne20]*rate[irnega] + y[mg24]*rate[irmgga]- y[ne20]*y[he4]*rate[irneag];
	
	dydt[mg24]=0.5*y[c12]*y[o16]*rate[ir1216]-y[mg24]*rate[irmgga]  + y[ne20]*y[he4]*rate[irneag]+y[si28]*rate[irsiga] - y[mg24]*y[he4]*rate[irmgag];
	
	dydt[si28]=0.5*y[c12]*y[o16]*rate[ir1216]+y[o16]*y[o16]*rate[ir1616] - y[si28]*rate[irsiga] + y[mg24]*y[he4]*rate[irmgag] - rate_irsi2ni*y[he4] + rate_irni2si*y[ni56];
	
	dydt[ni56]=rate_irsi2ni*y[he4] - rate_irni2si*y[ni56];
	}


/*
int main(){
	double y[N]={0.1,0.4,0.5,0,0,0,0};
	double rate[NRATE];
	double dfdy[N][N];
	double dydt[N];
	double temp=1e9;
	double dens=1e6;
	clock_t t=clock();
//	for(int i=0;i<1e6;i++){
		get_rate(temp,dens,y,rate);
		ydot(rate,y,dydt);
//		jacob(rate,y,dfdy);
//	}
//	t=clock()-t;
//	printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	for(int i=0;i<N;i++)
		cout<<dydt[i]<<endl;
	return(0);
}
*/	

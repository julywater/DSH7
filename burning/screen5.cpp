#include<iostream>
#include<math.h>
#include<time.h>
#include<stdio.h>
#define x13 (1.0/3.0) 
#define x14 (1.0/4.0)
#define x53 (5.0/3.0)
#define x532 (5.0/32.0)
#define x512 (5.0/12.0)
#define fact  1.25992104989487
#define co2  (x13 * 4.248719e3)
#define gamefx  0.3
#define gamefs  0.8
#define blend_frac  0.05
#define pi 3.1415926
#define SIZE 9
using namespace std;
void screen5(const double temp,const double den,const double zbar,const double abar,const double z2bar,double z1 [],double a1 [],double z2[],double a2[], double scor[],double scordt[],double scordd[])
{
	double zs13[SIZE],zs13inv[SIZE],zhat[SIZE],zhat2[SIZE],lzav[SIZE],aznut[SIZE];
	for(int jscreen=0;jscreen<SIZE;jscreen++){
	zs13[jscreen]=pow((z1[jscreen]+z2[jscreen]),x13);
	zs13inv[jscreen]=1.0/zs13[jscreen];
    zhat[jscreen]=pow((z1[jscreen]+z2[jscreen]),x53)-pow(z1[jscreen],x53)-pow(z2[jscreen],x53);
    zhat2[jscreen]= pow((z1[jscreen]+z2[jscreen]),x512)- pow(z1[jscreen],x512)-pow(z2[jscreen],x512);
    lzav[jscreen]= x53 * log(z1[jscreen]*z2[jscreen]/(z1[jscreen] + z2[jscreen]));
    aznut[jscreen]= pow((z1[jscreen]*z1[jscreen]*z2[jscreen]*z2[jscreen]*a1[jscreen]*a2[jscreen]/(a1[jscreen]+a2[jscreen])),x13);
	}
	//abar zaber dependence
	const double ytot=1.0/abar;
	const double rr=den*ytot;
	const double tempi=1.0/temp;
	const double dtempi=-tempi*tempi;
	const double deni =1.0/den;
	
	const double pp= sqrt(rr*tempi*(z2bar + zbar));
    double qq= 0.5/pp *(z2bar + zbar);
    const double dppdt= qq*rr*dtempi;
    const double dppdd=qq*ytot*tempi;

    const double qlam0z   = 1.88e8 * tempi * pp;
    const double qlam0zdt = 1.88e8 * (dtempi*pp + tempi*dppdt);
    const double qlam0zdd = 1.88e8 * tempi * dppdd;

    const double taufac   = co2 *pow(tempi,x13);
    const double taufacdt = -x13*taufac*tempi;

    qq=rr*zbar;
    const double xni=pow(qq,x13);
    const double dxnidd=x13 * xni * deni;

    const double aa= 2.27493e5*tempi*xni;
    const double daadt= 2.27493e5 * dtempi * xni;
    const double daadd= 2.27493e5 * tempi * dxnidd;
	//
	
	//individual screen
	for(int jscreen=0;jscreen<SIZE;jscreen++){
	double bb =z1[jscreen]* z2[jscreen];
    double gamp= aa;
    double gampdt= daadt;
    double gampdd= daadd;

    qq= fact * bb * zs13inv[jscreen];
    double gamef= qq * gamp;
    double gamefdt= qq * gampdt;
    double gamefdd  = qq * gampdd;

    double tau12 = taufac * aznut[jscreen];
    double tau12dt = taufacdt * aznut[jscreen];

    qq= 1.0/tau12;
    double alph12=gamef * qq;
    double alph12dt=(gamefdt - alph12*tau12dt) * qq;
    double alph12dd = gamefdd * qq;
	
//limit alph12 to 1.6 to prevent unphysical behavior.
//this should really be replaced by a pycnonuclear reaction rate formula
      if (alph12>1.6) {
       alph12   = 1.6;
       alph12dt = 0.0;
       alph12dd = 0.0;

       gamef    = 1.6 * tau12;
       gamefdt  = 1.6 * tau12dt;
       gamefdd  = 0.0;

       qq       = zs13[jscreen]/(fact * bb);
       gamp     = gamef * qq;
       gampdt   = gamefdt * qq;
       gampdd   = 0.0;
      }

	
//weak screening regime
      double h12w    = bb * qlam0z;
      double dh12wdt = bb * qlam0zdt;
      double dh12wdd = bb * qlam0zdd;

      double h12     = h12w;
      double dh12dt  = dh12wdt;
      double dh12dd  = dh12wdd;
	  
//intermediate and strong sceening regime
      if (gamef>gamefx){
		double gamp14   = pow(gamp,x14);
        double rr= 1.0/gamp;
        double qq= 0.25*gamp14*rr;
        double gamp14dt = qq * gampdt;
        double gamp14dd = qq * gampdd;
	  
		double cc = 0.896434* gamp * zhat[jscreen]\
                  - 3.44740* gamp14 * zhat2[jscreen]\
                  - 0.5551*(log(gamp) + lzav[jscreen])\
                  - 2.996;

        double dccdt = 0.896434* gampdt * zhat[jscreen]\
                  - 3.44740*gamp14dt * zhat2[jscreen]\
                  - 0.5551*rr*gampdt;

        double dccdd=0.896434* gampdd * zhat[jscreen]\
                  - 3.44740* gamp14dd * zhat2[jscreen]\
                  - 0.5551*rr*gampdd;
	
		
		double a3     = alph12 * alph12 * alph12;
		double da3    = 3.0* alph12 * alph12;

        qq     = 0.014+ 0.0128*alph12;
        double dqqdt  = 0.0128*alph12dt;
        double dqqdd  = 0.0128*alph12dd;

        rr     = x532 - alph12*qq;
        double drrdt  = -(alph12dt*qq + alph12*dqqdt);
        double drrdd  = -(alph12dd*qq + alph12*dqqdd);

        double ss     = tau12*rr;
        double dssdt  = tau12dt*rr + tau12*drrdt;
        double dssdd  = tau12*drrdd;

        double tt     =  -0.0098 + 0.0048*alph12;
        double dttdt  = 0.0048*alph12dt;
        double dttdd  = 0.0048*alph12dd;

        double uu     =  0.0055 + alph12*tt;
        double duudt  = alph12dt*tt + alph12*dttdt;
        double duudd  = alph12dd*tt + alph12*dttdd;

        double vv   = gamef * alph12 * uu;
        double dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt;
        double dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd;

        h12     = cc - a3 * (ss + vv);
        rr      = da3 * (ss + vv);
        dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt);
        dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd);

        rr     =  1.0- 0.0562*a3;
        ss     =  -0.0562*da3;
        drrdt  = ss*alph12dt;
        drrdd  = ss*alph12dd;
	     
		double xlgfac,dxlgfacdd,dxlgfacdt;
        if (rr >=0.77){
			xlgfac    = rr;
			dxlgfacdt = drrdt;
			dxlgfacdd = drrdd;}
        else{
			xlgfac    = 0.77;
			dxlgfacdt = 0.0;
			dxlgfacdd = 0.0;
        }
		
		h12    = log(xlgfac) + h12;
        rr     = 1.0/xlgfac;
        dh12dt = rr*dxlgfacdt + dh12dt;
        dh12dd = rr*dxlgfacdd + dh12dd;
		
		if (gamef<gamefs){
			rr     =  2.0*(gamefs - gamef);
			drrdt  = -2.0*gamefdt;
			drrdd  = -2.0*gamefdd;
			ss     = 2.0*(gamef - gamefx);
			dssdt  = 2.0*gamefdt;
			dssdd  = 2.0*gamefdd;
			
//store current values for possible blending
			double h12x    = h12;
			double dh12xdt = dh12dt;
			double dh12xdd = dh12dd;
			vv     = h12;
			h12    = h12w*rr + vv*ss;
			dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt;
			dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd;
			
// blend the transition region - from bill paxton
            if ((gamefs-gamef)<blend_frac*(gamefs - gamefx)){
				double alfa=(gamefs - gamef)/(blend_frac*(gamefs - gamefx));
				alfa   = 0.5 * (1.0- cos(pi*alfa));
				double beta= 1.0- alfa;
				h12= alfa * h12 + beta * h12x;
				dh12dt = alfa * dh12dt + beta * dh12xdt;
				dh12dd = alfa * dh12dd + beta * dh12xdd;
			}
       }
    }
	
//end of intermediate and strong screening if	
//machine limit the output
      h12    = max(min(h12,300.0),0.0);
      scor[jscreen]= exp(h12);
      if (h12 == 300.00){
       scordt[jscreen]= 0.0;
       scordd[jscreen]= 0.0;
	   }
      else {
       scordt[jscreen]= scor[jscreen]* dh12dt;
       scordd[jscreen]= scor[jscreen]* dh12dd;
      }
	  }
}



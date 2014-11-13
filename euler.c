
#include "../paul.h"
#include "../helm.h"
double get_om( double );
double get_om1( double );





static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static double explicit_viscosity = 0.0;
static int include_viscosity = 0;



void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   explicit_viscosity = theDomain->theParList.viscosity;
   include_viscosity = theDomain->theParList.visc_flag;
}

double get_omega( double * prim ){
   return( prim[UPP] );
}

void planetaryForce( struct planet * , int , double , double , double * , double * );
/*
void prim2tools( double * prim , double * Q , double r , double phi ){

   double rho = prim[RHO];
   Q[0] = rho;
   Q[1] = prim[PPP];
   double vr = prim[URR];
   double omega = prim[UPP];

   Q[2] = vr*r*(omega-pow(r,-1.5));
   Q[3] = rho*vr;

   double n = PHI_ORDER;
   double eps = G_EPS;
   double d  = sqrt( 1.+r*r-2.*r*cos(phi) );
   double dp = sin(phi);
   double fg = ( pow(d,n-1.)/pow( pow(d,n) + pow(eps,n) ,1.+1./n) );

   Q[4] = (rho-1.)*r*r*sin(phi)/pow( 1.+r*r-2.*r*cos(phi) + G_EPS*G_EPS , 1.5);
   Q[5] = (rho-1.)*r*(fg*dp/d)*(2.*M_PI*r);

   double rH = pow(PLANET_MASS/3.,1./3.)/DISK_MACH;
   
   double p = 0.8;
   double fd = 1./(1.+exp(-( d/rH - p )/(p/10.) ) );
   Q[6] = Q[5]*fd;
   
}*/

void prim2cons( double * prim , double temp,double * cons , double r , double dV ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double om  = get_om( r );
   double vp_off = vp - om*r;

   double v2  = vr*vr + vp_off*vp_off + vz*vz;
   double rhoe;
   double cs;
   invert_helmeos_pd(Pp,rho,&temp,&rhoe,&cs);
//   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe )*dV;
   cons[SRR] = rho*vr*dV;
   cons[LLL] = r*rho*vp*dV;
   cons[SZZ] = rho*vz*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( double * prim ,double temp, double * Ustar , double r , double Sk , double Ss , double * n , double * Bpack ){

   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];

   double om = get_om( r );
   double vp_off = vp - om*r;
   double v2 = vr*vr+vp_off*vp_off+vz*vz;

   double vn = vr*n[0] + vp*n[1] + vz*n[2];

   double vn_off = vn - om*r*n[1];
   double Ss_off = Ss + vn_off - vn;

 //  double rhoe = Pp/(gamma_law - 1.);
   double rhoe=0.0;
   double cs=0.0;
   invert_helmeos_pd(Pp,rho,&temp,&rhoe,&cs);



   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] =   rhostar*( vr + (Ss-vn)*n[0] );
   Ustar[LLL] = r*rhostar*( vp + (Ss-vn)*n[1] );
   Ustar[SZZ] =   rhostar*( vz + (Ss-vn)*n[2] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( double * cons , double * prim , double *temp ,double r , double dV ){
   
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/dV;
   double Sp  = cons[LLL]/dV/r;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( r );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om*r;
   double vz = Sz/rho;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + Sz*vz );
   double rhoe = E-KE;
   
   double Pp=0;
   double cs=0;
   invert_helmeos_ed(rhoe,rho,temp,&Pp,&cs);

//   printf("%f   %f   %f\n",*temp,rhoe,Pp);

//   double Pp = (gamma_law - 1.)*rhoe;

   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;




   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = vp/r;
   prim[UZZ] = vz;
  
//printf("%f   %f   %f\n",*temp,rhoe,Pp);
 
   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}


	

void flux( double * prim , double * flux , double r , double * n ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double om  = get_om( r );

   double vn = vr*n[0] + vp*n[1] + vz*n[2];
   double wn = om*r*n[1];
   double vp_off = vp - om*r;

   double rhoe = Pp/(gamma_law - 1.);
   double v2 = vr*vr + vp_off*vp_off + vz*vz;

   flux[DDD] = rho*vn;
   flux[SRR] = rho*vr*vn + Pp*n[0];
   flux[LLL] = r*rho*vp*vn + r*Pp*n[1];
   flux[SZZ] = rho*vz*vn + Pp*n[2];
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*wn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

double get_dp( double , double );

void source( double * prim , double * cons , double * xp , double * xm , double dVdt ){
   
   double rp = xp[0];
   double rm = xm[0];
   double dphi = get_dp(xp[1],xm[1]);
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double r_1  = .5*(rp+rm);
   double r2_3 = (rp*rp + rp*rm + rm*rm)/3.;
   double vr  = prim[URR];
   double omega = prim[UPP];
 
   double centrifugal = rho*omega*omega*r2_3/r_1*sin(.5*dphi)/(.5*dphi);
   double press_bal   = Pp/r_1;

   cons[SRR] += dVdt*( centrifugal + press_bal );

   double om  = get_om( r_1 );
   double om1 = get_om1( r_1 );

   cons[TAU] += dVdt*rho*vr*( om*om*r2_3/r_1 - om1*(omega-om)*r2_3 );
 
   if( include_viscosity ){
      double nu = explicit_viscosity;
      cons[SRR] += -dVdt*nu*rho*vr/(r_1*r_1);
   }

}

void visc_flux( double * prim , double * gprim , double * flux , double r , double * n ){

   double nu = explicit_viscosity;

   double rho = prim[RHO];
   double vr  = prim[URR];
   double om  = prim[UPP];
   double om_off = om - get_om(r);
   double vz  = prim[UZZ];

   double dnvr = gprim[URR];
   double dnom = gprim[UPP];
   double dnvz = gprim[UZZ];

   flux[SRR] = -nu*rho*( dnvr - n[1]*2.*om );
   flux[LLL] = -nu*rho*( r*r*dnom + n[1]*2.*vr );
   flux[SZZ] = -nu*rho*dnvz;
   flux[TAU] = -nu*rho*( vr*dnvr+r*r*om_off*dnom+vz*dnvz );//- 2.*r*om_off*om );

}

void vel( double * prim1 , double temp1,double * prim2 ,double temp2, double * Sl , double * Sr , double * Ss , double * n , double r  ){
	
   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r + prim1[UZZ]*n[2];
   double cs1=0.0;
 //  double temp1=1.25e9;
   double rhoe1=0.0;
   invert_helmeos_pd(P1,rho1,&temp1,&rhoe1,&cs1);


//   double cs1 = sqrt(gamma_law*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r + prim2[UZZ]*n[2];
   
   double cs2=0.0;
//   double temp2=1.25e9;
   double rhoe2=0.0;
   invert_helmeos_pd(P2,rho2,&temp2,&rhoe2,&cs2);

//   double cs2 = sqrt(gamma_law*P2/rho2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}

double get_dL( double * , double * , int );

double mindt(double * prim , double temp, double w , double * xp , double * xm ){

   double r = .5*(xp[0]+xm[0]);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r;
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   
   double eout,cs;
   invert_helmeos_pd(Pp,rho,&temp,&eout,&cs);

   double maxvr = cs + fabs(vr);
   double maxvp = cs + fabs(vp);
   double maxvz = cs + fabs(vz);

   double dtr = get_dL(xp,xm,0)/maxvr;
   double dtp = get_dL(xp,xm,1)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;

   return( dt );

}

double getReynolds( double * prim , double w , double r , double dx ){

   double nu = explicit_viscosity;

   double vr = prim[URR];
   double omega = prim[UPP];
   double vp = omega*r-w;
   double vz = prim[UZZ];

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double cs = sqrt(gamma_law*Pp/rho);
   
   double v = sqrt(vr*vr + vp*vp + vz*vz);

   double Re = (v+cs)*dx/nu;
   
   return(Re);

}


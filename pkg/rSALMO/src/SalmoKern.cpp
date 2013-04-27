/*****************************************************************************/
/* Kernel module of SALMO 
   contains separate functions for advection, diffision and reaction
   call be called from R or any other integration interface

   JAVA version:         Steffen Dietze TU Ilmenau, 1997 
                         with bugfixes and extensions of thpe
   C version             S. Rolinski (SR)      2004 - 2008  
                         T. Petzoldt (thpe)    2008 -

   Thanks to: René Sachse (and many more, who just need to be asked if they
              want to be mentioned here)
*/
/*****************************************************************************/


/* uncomment the following line if you want to call this .DLL/.so from R     */
/* and want to print to screen with Rprintf                                  */
//#define use_R true

#ifdef use_R
#include <R.h>                    /* for R compatibility */
#include <Rinternals.h>           /* for R compatibility */
#endif

#ifndef math_h
#include <math.h>                 /* needed when R headers are not used     */
#endif

#include <algorithm>
using namespace std;

//#include <iostream.h>           /* only necessary for specific test cases  */
#include <iostream> 
#include <stdio.h>                /* needed for debugging with printf        */

/*===========================================================================*/
/*  General Settings, Conventions and Helper Formulae                        */
/*===========================================================================*/

/* Positions of Constants and Parameters in Vectors c, p und u               */
#include "salmo_constants.h" 

/* Macros to enhance readability of dynamic data structures                  */
/*  (parameters and temporary variables for phytoplankton)                   */
#include "salmo_macros.h"

#include "salmo_headers.h"


double saturation(double T) {
  if (T < 0) T = 0;
  T += 273.15;  // Celsius to Kelvin
  return exp(-139.34411+157570.1/T-66423080./
    pow(T,2)+12438000000./pow(T,3)-862194900000./pow(T,4));
}


/******************************************************************************/
/*  S A L M O  Kernel                                                         */
/*  can be called by 'all_layers'                                             */
/******************************************************************************/


void SalmoKern(int* nOfVar, double* c, double* p, double* u, double* x, double* dxq, double* dxs) {
 int np = nOfVar[inumberOfParameters];    // important for accessing "p"
 int nx = nOfVar[inumberOfPhyto];         // number of phytoplankton groups: 3
 
 /* derivatives of the state variables, splitted in source and sink terms     */
 double dNq = 0.0;
 double dNs = 0.0;
 double dPq = 0.0;
 double dPs = 0.0;
 double dZq = 0.0;
 double dZs = 0.0;
 double dDq = 0.0;
 double dDs = 0.0;
 double dOq = 0.0;
 double dOs = 0.0;
 
 /* phytoplankton and carbon:chl ratio G can have variable number nx          */
 double *X   = new double[nx];
 double *dXq = new double[nx];
 double *dXs = new double[nx];
 double *G   = new double[nx];
 double *dG  = new double[nx];
 double *xmineralisation = new double[nx];
 int j, zj;

 /* state variables in the following order                                    */
 double N   = x[0];
 double P   = x[1];
 for (j = 0; j < nx; j++) X[j]  = x[2 + j];
 double Z   = x[2 + nx];
 double D   = x[3 + nx];
 double O   = x[4 + nx];
 for (j = 0; j < nx; j++) G[j]  = x[4 + nx + j];

 /*===========================================================================*/
 /* Assignment of parameters from vectors 'p' and 'c'                         */
 /*                                                                           */ 
 /* Notes:                                                                    */
 /*   - the cXXXXX constants are defined in salmo_constants.h                 */
 /*   - integration time is also contained in vector 'c'                      */
 /*   - all the XXX_j are #define macros pointing to the dynamic array pC     */
 /*   - loops over XXX_j must have loop variable 'j' !!!                      */
 /*   - see salmo_macros.h                                                    */
 /*===========================================================================*/

 double EPSMIN  = c[cEPSMIN];
 double ZLIGHT  = c[cZLIGHT];
 int CHLC       = (int) c[cCHLC];
 double Zres    = c[cZres];    
 double SF      = c[cSF];
 int simyear    = (int) c[csimyear];
 double ANSFMIN = c[cANSFMIN]; 
 double APSFMAX = c[cAPSFMAX]; 
 double APSFMIN = c[cAPSFMIN];
 double APSFT   = c[cAPSFT];   
 double AZMAX   = c[cAZMAX];   
 double AZMIN   = c[cAZMIN];
 double DTA     = c[cDTA];     
 double DTB     = c[cDTB];     
 double DTC     = c[cDTC];
 double DTMIN   = c[cDTMIN];
 double EPSD    = c[cEPSD];  
 //double GI      = c[cGI];
 double GMAX    = c[cGMAX];    
 double GMIN    = c[cGMIN];
 double KANSF   = c[cKANSF];   
 double KAPSF   = c[cKAPSF];    
 double KDEN    = c[cKDEN];
 double KNDST   = c[cKNDST];   
 double KMINER  = c[cKMINER];
 double KMO     = c[cKMO];     
 //double KPSED   = c[cKPSED];    
 double KSEZA   = c[cKSEZA]; 
 double KSRFMAX = c[cKSRFMAX]; 
 double KXG     = c[cKXG];      
 double KXMIN   = c[cKXMIN]; 
 double KZMIN   = c[cKZMIN];
 double LGH     = c[cLGH];     
 double LGL     = c[cLGL];      
 double LINDEN  = c[cLINDEN];
 double LXH     = c[cLXH];     
 double LXHN    = c[cLXHN];
 double LXL     = c[cLXL];     
 double LXLN    = c[cLXLN];
 double MGH     = c[cMGH];     
 double MGL     = c[cMGL];
 double MOMIN   = c[cMOMIN];   
 double MOT     = c[cMOT];
 double MXH     = c[cMXH];     
 double MXL     = c[cMXL];
 double NDSSTART = c[cNDSSTART]; 
 double NDSEND = c[cNDSEND];  
 double NDSMAX = c[cNDSMAX];
 double KNDS    = c[cKNDS];    
 double OPTNP   = c[cOPTNP];
 double PF      = c[cPF];
 double R       = c[cR];
 double RAT     = c[cRAT];     
 double RATF    = c[cRATF];
 double RATN    = c[cRATN];    
 double RATNF   = c[cRATNF];
 //double RL      = c[cRL];    
 //double RLW     = c[cRLW];
 double RXMF    = c[cRXMF]; 
 double RZMIN   = c[cRZMIN];   
 double RZOPT   = c[cRZOPT];    
 double RZTMIN  = c[cRZTMIN];
 double EPSR    = c[cEPSR];
 double TOPTZ   = c[cTOPTZ];
 double UXZD    = c[cUXZD];
 double VD      = c[cVD];    
 //double VMIG    = c[cVMIG];
 double WPKX    = c[cWPKX];    
 double WPKZ    = c[cWPKZ];
 double YD      = c[cYD];      
 double YND     = c[cYND];      
 double YNX     = c[cYNX]; 
 double YOX     = c[cYOX];     
 double YZN     = c[cYZN];
 double YZP     = c[cYZP];
 // bool densed    = false;
 double SEZMAX  = c[cSEZMAX];
 double KO      = c[cKO];
 // double VO2ATM = c[cVO2ATM];      // atmospaeric re-aereation (mg/m)
 int npsfmode  = c[cNpsfmode];    // true: P release only nitrate dependent (shallow lakes)
 int srfmode   = c[cSrfmode];     // true: srf column is strong rain factor; false: is epsmin
 int xminermode = c[cXminermode]; // way to calculate oxygen depletion of X and D sedimentation


// --- testing code ---
// printf("- (%d) (%e)", xminermode, VO2ATM);

 double *pC       = new double[np*nx];
 int    phyn;
 for (j = 0; j < nx; j++) {
   phyn       = (int) p[j*np];
   pC[j]      = p[j*np];               //?? ThPe, check this
   EPSX_j     = p[phyn*np + pEPSX];
   BETA_j     = p[phyn*np + pBETA];
   GAMMAX_j   = p[phyn*np + pGAMMAX];
   GAMMIN_j   = p[phyn*np + pGAMMIN];
   KI_j       = p[phyn*np + pKI];
   KN_j       = p[phyn*np + pKN];
   KP_j       = p[phyn*np + pKP];
   KPF_j      = p[phyn*np + pKPF];
   NFIX_j     = p[phyn*np + pNFIX];
   PFC_j      = p[phyn*np + pPFC];
   PFX_j      = p[phyn*np + pPFX];
   PHOTXMAX_j = p[phyn*np + pPHOTXMAX];
   PHOTXMIN_j = p[phyn*np + pPHOTXMIN];
   RXTMIN_j   = p[phyn*np + pRXTMIN];
   RXTOPT_j   = p[phyn*np + pRXTOPT];
   RXTSLOPE_j = p[phyn*np + pRXTSLOPE];
   TOPTX_j    = p[phyn*np + pTOPTX];
   UXZ_j      = p[phyn*np + pUXZ];
   VS_j       = p[phyn*np + pVS];
   YX_j       = p[phyn*np + pYX];
   ZETA_j     = p[phyn*np + pZETA];
 }
 

 /* Dynamic Vector for intermediate values for each phytoplankton group
    Note assignments in salmo_macros.h                                        */
 double *cB  = new double[nx*(cB_rows+1)];

 /*===========================================================================*/
 /* Assign elements of input matrix                                           */
 /* here "u" contains data for one layer with depth information               */
 /* In R use:                                                                 */
 /*   u1 <- cbind(t, v, zmixreal, zmix, qin, qout,                            */
 /*     srf, iin, te, nin, pin, pomin, x1in, x2in, zin, oin, ae, ah)          */
 /*===========================================================================*/

 double t        = u[undt];
 double v        = u[uvol];
 double depth    = u[udepth]; // water depth of lower boundary of layer
 double zmix     = u[uzmix];  // layer depth
 double qin      = u[uqin];
 double ased     = u[uased];  // former qout is now sediment contact area of layer
 double srf      = u[usrf];
 double iin      = u[uiin];
 double temp     = u[utemp];
 double nin      = u[unin];
 double pin      = u[upin];
 double pomin    = u[upomin];
 for (j = 0;j < nx;j++) { cB[j+ixin*nx] = u[ux1in + j]+X[j]; }
 double zin      = u[uzin];
 double oin      = u[uoin];
 double aver     = u[uaver];

 /*===========================================================================*/
 /* Helper equations for ODE system                                           */
 /*===========================================================================*/  

 /*---------------------------------------------------------------------------*/
 /* 15.0    Oxygen during full circulation and in epilimnion                  */
 /*       = saturation concentration                                          */
 /*---------------------------------------------------------------------------*/
 // double sat = saturation(temp);                                    /* 15.1 */ 
 // double o = sat; // thpe: redundant in the recent version

 /*---------------------------------------------------------------------------*/
 /* Calculate nkons                                                           */
 /*---------------------------------------------------------------------------*/

 int M   = (int) (zmix / ZLIGHT);                                      /* 9.3 */
 if (M < 2) M = 2;
 
 double eps,ksrf,kxn,kx,iredz,nkons,nresp;
 double photxj = 0;   // ToDo: rename this to avoid confusion with photx_j
 double zvonj  = 0;
 // RL now considered externally
 // double ired   = RL * iin; RL: Reflection at surface                /* 7.1 */
 double ired   = iin;                                                  /* 7.1 */

 for (j = 0;j < nx;j++) {
   wx_j = 0;  rx_j = 0;  photx_j = 0.0;  dG[j] = 0.0;
 }

 /* !srfmode means that measured epsilon is provided in srf data column       */
 if (!srfmode) {
   eps = srf;
 } else {
   ksrf = (1.5 + 0.5 * cos((t - 30) * PI * 2 / simyear)) * srf;       /* 7.4  */  
   if ((qin > 0) && (ksrf > KSRFMAX)) {                               /* 7.3  */
     eps = EPSMIN + EPSR * (ksrf - KSRFMAX); 
   } else { 
     eps = EPSMIN; 
   }
 }  

 /* light extinction due to phytoplankton and detritus                        */
 double epsxi = 0;
 for (j = 0;j < nx;j++) epsxi = epsxi + EPSX_j * X[j];
 eps = (eps + epsxi + EPSD * D);                                      /* 7.2  */ 

 for (j = 0;j < nx;j++) {
   phoxt_j = (PHOTXMAX_j-PHOTXMIN_j) / TOPTX_j * temp + PHOTXMIN_j;   /* 9.7  */
 }

 if (N > (WPKX * OPTNP))                                              /* 9.12 */
   kxn = KXMIN + LXHN * pow(N, MXH);
 else
   kxn = LXLN * pow(N, MXL); 

 /* if(kxn == 0)    Protection against division by zero */            /* 9.11 */
 /* phoxn_j  = N / X[j] / ((KN_j + N) / X[j]);
 else   */
 for (j = 0; j < nx; j++)
   phoxn_j = N / X[j] / ((KN_j + N) / kxn + (KN_j + N) / X[j]);       /* 9.11 */

 if (P > WPKX)                                                        /* 9.10 */
   kx = KXMIN + LXH * pow(P, MXH); 
 else
   kx = LXL * pow(P, MXL); 


 /* if(kx == 0)     Protection against division by zero */            /* 9.9  */
 /* phoxp_j  = P / X[j] / ((KP_j+P) / X[j]); 
 else  */
 for (j = 0; j < nx; j++)
   phoxp_j    = P / X[j] / ((KP_j + P) / kx + (KP_j + P) / X[j]);     /* 9.9  */


 /* N-Limitation or P-Limitation or N Fixation                                */
 for (j = 0;j < nx;j++) {
   if ((N / P < OPTNP) && (NFIX_j < 0.0001))                          /* 9.6  */
     phoxns_j = phoxn_j;
   else 
     phoxns_j = phoxp_j;
   rxt_j = (RXTOPT_j - RXTMIN_j) / TOPTX_j * temp + RXTMIN_j;         /* 9.14 */
 }

 /*---------------------------------------------------------------------------*/
 /* Primary production according to the traditional formula (daily average)   */
 /* or Baumert's Chl:C ratio method                                           */
 /* Note: the Chl:C method is not yet completely implemented                  */
 /*---------------------------------------------------------------------------*/
 
 if (! CHLC) {
   for (zj = 1; zj <= M; zj++) {                                      /* 9.2  */
     zvonj = zmix * (zj - 1) / (double) M;                            /* 7.5  */
#ifdef use_R
     Rprintf("zj = %d, M = %d, zvonj = %8.4f, zmix = %8.2f\n", zj, M, zvonj, zmix);
#endif
     /* 7.0 Light in the mixed layer */
     iredz = ired * exp(-eps * zvonj);                                /* 7.0  */
     for (j = 0; j < nx; j++) {
       if ((N / P) < OPTNP)                                           /* 9.8  */
         phoxl_j= iredz * (1-NFIX_j) / (KI_j + iredz);
       else 
         phoxl_j= iredz / (KI_j + iredz);                             /* 9.8  */

       // !!! thpe: photxj vs. photx_j is confusing here; rename variables
       photxj = phoxt_j * phoxl_j * phoxns_j;                         /* 9.5  */
       photx_j = photx_j + photxj;
       // rx_j  = rxt_j + RXMF * photx_j;                             /* 9.13 */
       // wx_j  = wx_j + photx_j - rx_j;
#ifdef use_R
       //Rprintf("%g\n", photx_j);
#endif
     } /* end j  (loop over all phytoplankton groups)                         */
   } /* end zj  (loop over M sub-layers)                                      */
 } else {
    double alpha, phox0, phoxlj, dGj;
    for (zj = 1; zj <= M; zj++) {                                     /* 9.2  */
      zvonj = (zj - 1) / M * zmix;                                    /* 7.5  */

      /* Light in the mixed layer */
      iredz = ired * exp(-eps * zvonj);                               /* 7.0  */

      for (j = 0; j < nx; j++) {
        phox0 = phoxt_j * phoxns_j;

        G[j] = 1.0; // thpe: still a workaround, sets CHL:C to one

        // thpe: check parameter 8.64 ??  day --> seconds
        alpha = phox0 / (KI_j / 8.64 * G[j]);
        phoxlj = phox0 * (1 - exp(-alpha * G[j] * iredz / phox0) *
                   exp(-BETA_j * alpha * G[j] * iredz / phox0));
        dGj   = G[j] * phox0 * (exp(-alpha * G[j] * iredz / phox0)-
                  (G[j] - GAMMIN_j) / (GAMMAX_j - GAMMIN_j)) *
                  exp(-BETA_j * alpha * G[j] * iredz / phox0);
                
        // photxj = phoxt_j * phoxl_j * phoxns_j;                      /* 9.5  */
        photx_j = photx_j + phoxlj;
        dG[j]   = dG[j] + dGj;       
      } // end j
    } //end zj (loop over M)
 }
 /* End of branching of photosynthesis sub-models */

 /*---------------------------------------------------------------------------*/
 /* Light at the bottom of the layer                                          */
 /*---------------------------------------------------------------------------*/
 double ibottom = ired * exp(-eps * zmix);  

 for (j = 0;j < nx;j++) { 
   photx_j = photx_j / M; 
   dG[j] = dG[j] / M; 
   rx_j =   rxt_j + RXMF * photx_j;                                   /* 9.13 */
   wx_j =    wx_j + photx_j - rx_j;
 }   

 nkons = 0;
 nresp = 0;
 if ((N / P) < OPTNP) {                                                /* 1.2 */
   for (j = 0;j < nx;j++) { 
     if (NFIX_j < 0.0001) { 
       nkons = nkons + photx_j * X[j]; 
       nresp = nresp + rx_j * X[j]; 
     } 
   }
 } else {
   for (j = 0;j < nx;j++) { 
     nkons = nkons + photx_j * X[j]; 
     nresp = nresp + rx_j * X[j]; 
   } 
 }
 nkons = nkons / YNX;
 nresp = nresp / YNX;
 
 /*---------------------------------------------------------------------------*/
 /* Calculations for nrem (remineralisation of nitrogen)                      */
 /*---------------------------------------------------------------------------*/

 double kzd, hwgdbd, hwgdd, sumxpf, kz, gdb, gd, hwgd, g, g3d, hwgi;
 double az, assiz, rzt, rzg, rz, wz, egg, nexkr, nrem, dt;

 double mortz = (MOMIN + MOT * temp) * Z / (KMO + Z);                /* 12.12 */
 double nmort = mortz * RATNF / YZN;                                 /*  1.4  */
 double gdt = (GMAX - GMIN) * exp(-R * fabs(log(temp / TOPTZ))) + GMIN; /* 9.22 */
 
 // ToDo: Reformulate remaining occurences of cB and pC 

 cB[nx-1 + nx*ipf] = pC[nx-1 + nx*pPFC];                              /* 9.25 */
 double xsum = 0;
 for (j = (nx-2); j >= 0; j--) {
   xsum = xsum + cB[j + 1 + nx*ipf] * X[j + 1];
   if(xsum > PFX_j)                                                   /* 9.25 */
     pf_j = 1/(xsum-PFX_j) / (1/(KPF_j-PFX_j) + 1/(xsum-PFX_j));
   else 
     pf_j = 1;
 }
 
 if((D * PF) > WPKZ)                                                  /* 9.24 */
   kzd = KZMIN + LGH * pow((D * PF), MGH); 
 else 
   kzd = LGL * pow((D * PF), MGL);

 for (j = 0;j < nx;j++) {
   if ((X[j] * pf_j) > WPKZ)                                           /* 9.24 */
     kzi_j = KZMIN + LGH * pow((X[j] * pf_j), MGH); 
   else 
     kzi_j = LGL * pow((X[j] * pf_j), MGL);
 }

 /* if(kzd == 0)  protects against division by zero */  
 /* hwgdbd= D * PF / Z / ((KXG + D) / Z);
 else */
 hwgdbd= D * PF / Z / ((KXG + D) / kzd + (KXG + D) / Z);             /* 9.23 */

 for (j = 0;j < nx;j++) {
   /* if(kzi_j == 0)     protects against difvision by zero */  
   /* hwgdb_j= X[j] * pf_j / Z / ((KXG + X[j]) / Z);
   else */
   hwgdb_j= X[j] * pf_j / Z / ((KXG + X[j]) / kzi_j + (KXG + X[j]) / Z); /* 9.23 */
 }

 hwgdd  = gdt * hwgdbd;                                               /* 9.21 */
 sumxpf = 0;
 for (j = 0;j < nx;j++) {
  hwgdi_j  = gdt * hwgdb_j;                                           /* 9.21 */
  sumxpf =  sumxpf + X[j] * pf_j;
 }
 sumxpf = sumxpf + D * PF; /* helper variable */

 if(sumxpf > WPKZ)                                                    /* 9.31 */
   kz  = KZMIN + LGH * pow(sumxpf, MGH);
 else 
   kz  = LGL * pow(sumxpf, MGL);

 /* if(kz == 0)  protects against division by zero */
 /* gdb = sumxpf / Z / ((KXG + sumxpf) / Z);
 else */
 gdb   = sumxpf / Z / ((KXG + sumxpf) / kz + (KXG + sumxpf) / Z);     /* 9.30 */
 gd    = gdt * gdb;                                                   /* 9.29 */
 hwgd  = hwgdd;                                                       /* 9.20 */
 for (j = 0; j < nx; j++) hwg_j = hwgdi_j;                            /* 9.20 */
 
 g    = gd;                                                           /* 9.28 */
 hwgi = 0;
 for (j = 0;j < nx;j++) hwgi = hwgi + hwg_j;

 g3d = hwgd * g / (hwgi + hwgd);                                      /* 9.19 */

 for (j = 0;j < nx;j++) gi_j = hwg_j * g / (hwgi + hwgd);             /* 9.19 */

 az      = AZMAX - (AZMAX - AZMIN) * g / GMAX;                        /* 12.6  */

 assiz   = 0;
 for (j = 0;j < nx;j++) assiz = assiz + gi_j * UXZ_j;

 assiz   = az * (assiz + g3d * UXZD);                                 /* 12.5  */
 rzt     = (RZOPT - RZTMIN) * pow(temp/ TOPTZ,2) + RZTMIN;            /* 12.9  */
 rzg     = ((RZOPT - RZMIN) * g / GMAX + RZMIN) / RZOPT;              /* 12.8  */
 rz      = rzg * rzt;                                                 /* 12.7  */
 wz      = assiz - rz;                                                /* 12.4  */
 dt      = exp(DTA -DTB * log(temp) -DTC * log(temp) * log(temp));    /* 12.3  */

 if (dt < DTMIN)                                                      /* 12.2  */
   egg   = 1;
 else 
   egg   = DTMIN / dt;

 wz      = egg * wz;
 assiz   = egg * assiz;
 rz      = egg * rz;

 nexkr = 0;
 for (j = 0; j < nx; j++) nexkr = nexkr + gi_j; 

 nexkr   = ( nexkr / YNX + g3d / YND + rz / YZN ) * RATN;             /* 1.5  */
 nkons   = nkons + assiz / YZN * RATN * Z;            
 nrem    = (nmort + nexkr) * Z;                                       /* 1.3  */

 /*---------------------------------------------------------------------------*/
 /* N release resp. denitrification at the sediment surface                   */
 /*---------------------------------------------------------------------------*/
 double ansfq,ansfs;
 if ((NDSSTART <= fmod(t,simyear)) && (fmod(t,simyear) < NDSEND)) {
   ansfs = NDSMAX  * N /(KNDS + N) * pow(KNDST, temp-4);
   ansfq = ANSFMIN;
 } else { 
   ansfq = ANSFMIN + KANSF * temp;                                    /* 1.8  */
   ansfs = 0.0; 
 } 

 /*---------------------------------------------------------------------------*/
 /* Calculations for pkons, prem, psed                                        */
 /*---------------------------------------------------------------------------*/
 double pkons, pmort, pexkr, presp, prem; /*,psed;*/
 pkons = 0;
 pexkr = 0;
 presp = 0;
 for (j = 0;j < nx;j++) {                                             /* 4.3  */
   //  pkons = pkons + wx_j * X[j] / YX_j;
   pkons = pkons + photx_j * X[j] / YX_j;        // Sink
   presp = presp + rx_j * X[j] / YX_j;           // Source
   pexkr = pexkr + gi_j / YX_j;                  // Source
 }
 
 pmort   = mortz * RATF / YZP;                                        /* 4.5  */

 // Separation in sources and sink terms
 // pexkr   = (pexkr + g3d / YD - wz / YZP) * RAT;                    /* 4.6  */
 pexkr   = (pexkr + g3d / YD + rz / YZP) * RAT;  // Source            /* 4.6  */
 pkons   = pkons + assiz / YZP * RAT * Z;  // Sink                    
 prem    = (pmort + pexkr) * Z;            // Source                  /* 4.4  */
 // psed    = KPSED * P;                                              /* 4.8  */


 /*---------------------------------------------------------------------------*/
 /* P-release dependend on O2 or NO3-                                         */
 /*---------------------------------------------------------------------------*/
 double apsf; // always positive, therefore no source/sink splitting


 if (npsfmode == false) { /* default: O2 AND NO3- dependent                   */
   double oplusn = O + N/0.3;
   if(oplusn <= LINDEN) {                                             /* 4.10 */
     apsf = APSFMAX; 
   } else {
     apsf = APSFMAX /(oplusn - LINDEN) /(1 /(KAPSF - LINDEN) + 1/(oplusn - LINDEN)) + APSFMIN; 
   }
 } else {/* only NO3-dependent                                                */
   double klinden = 0.3 * LINDEN;
   if(N <= klinden) {                                                 /* 4.10 */
     apsf = APSFMAX; 
   } else { 
     apsf = APSFMAX /(N-klinden) /(1 /(KAPSF - klinden) + 1/(N - klinden)) + APSFMIN; 
   }
 } 
 
 // Rprintf("APSFMIN, APSFMAX, apsf = %g %g %g\n", APSFMIN, APSFMAX, apsf);

 /*---------------------------------------------------------------------------*/
 /*  Temperature correction. Reference temperature = 4 deg C                  */
 /*---------------------------------------------------------------------------*/
 apsf = apsf * pow(APSFT, temp - 4);

 /*---------------------------------------------------------------------------*/
 /* Calculations for xwa, xsed, xgraz                                         */
 /*---------------------------------------------------------------------------*/
 for (j = 0;j < nx;j++) {
   xwa_j = wx_j * X[j];                                               /* 9.1  */

   // bx[j]    = VS[j] / zmix;

   /* Zres: resuspension depth */
   if (depth > Zres) 
     bx_j = VS_j/zmix * SF * (1 - aver);                              /* 9.17 */
   else 
     bx_j = 0;
  
   xsed_j  = bx_j * X[j];                                             /* 9.16 */
   xgraz_j = gi_j * Z;                                                /* 9.18 */
 }

 /*---------------------------------------------------------------------------*/
 /* Calculations for zwa, zmo                                                 */
 /*---------------------------------------------------------------------------*/
 // double zwa    = wz * Z;                                          /* 12.1  */
 double zmo    = mortz * Z;                                          /* 12.11 */
 // double miraz  = VMIG / zmix;                                     /* 13.2  */
 // double zmig   = miraz * Z;                                       /* 13.1  */


 /*---------------------------------------------------------------------------*/
 /* Calculations for dsed, dgraz                                              */
 /*---------------------------------------------------------------------------*/
 double bd       = 0;
 if (depth > Zres) 
   bd     = VD / zmix * SF * (1 - aver);                             /* 18.3  */

 double dsed   = bd * D;                                             /* 18.2  */
 double dgraz  = g3d * Z;                                            /* 18.4  */

 /*---------------------------------------------------------------------------*/
 /* Calculations for nim, nex, nsf,                                           */
 /*         pim, pex, psf,                                                    */
 /*         xex, xim,                                                         */
 /*         zim,                                                              */
 /*         dim, dex,                                                         */
 /*---------------------------------------------------------------------------*/
 double nim, nsfq, nsfs, pim, psf, zim, dim, dilution;

 dilution = qin / (v + qin);

 /* This formulation of import takes only the source,                         */
 /* while the sinks are directly in the final equation. SR 2003-12-03         */
 nim   = nin * dilution;                                              /* 1.1  */
 nsfq  = ansfq * ased / v;                                            /* 1.7  */
 nsfs  = ansfs * ased / v;                                            /* 1.7  */
 pim   = pin * dilution;                                              /* 4.1  */
 psf   = apsf * ased / v;                                             /* 4.9  */

 for (j = 0;j < nx;j++) xim_j = xin_j * dilution;                     /* 9.34 */
 
 zim  = zin * dilution;                                              /* 13.13 */
 dim  = pomin * dilution;                                            /* 18.1  */


 /*---------------------------------------------------------------------------*/
 /* SR 16.06.2003: former mixing calculation using "fluk"                     */
 /* is now found in the previously called function "advection"                */
 /* and then added to the results vector dx                                   */
 /*---------------------------------------------------------------------------*/
 // nfluk   = (NH - N) * ah / ve;                                     /*  2.2 */
 // pfluk   = (PH - P) * ah / ve;                                     /*  5.2 */
 // for (j = 0;j < nx;j++) { xfluk[j]  = (XH[j] - X[j]) * ah / v; }   /* 10.1 */
 // zfluk   = (ZH - Z) * ah / ve;                                     /* 13.3 */
 // dfluk   = (DH - D) * ah / ve;                                     /* 19.2 */

 /*---------------------------------------------------------------------------*/
 /*  Calculations for oim, oex, oprod, okons                                  */
 /*---------------------------------------------------------------------------*/
 double oim, oprod, oresp, seza, lsez, minerd, sumxsed, lo, okons, olim;
 
 /* limitation of respiration by oxygen concentration (Rene Sachse 2012-07-09 */
 olim = O / ( KO + O );

// oim = (oin - O) * qin / (v + qin);                                /* 17.1  */
 oim = oin * dilution;                                               /* 17.1  */

// oprod = 0;
// for (j = 0;j < nx;j++) { oprod   = oprod + wx_j * X[j]; }
// oprod = oprod / YOX;                                              /* 17.3  */

 oprod = 0;
 oresp = 0;
 for (j = 0;j < nx;j++) { 
   oprod   = oprod + photx_j * X[j]; 
   oresp   = oresp + rx_j * X[j] * olim; 
 }
 oprod = oprod / YOX;                                                /* 17.3  */
 oresp = oresp / YOX;                                                /* 17.3  */

 // hidden Parameter 0.2 - 0.4,  SR 8.2.2005
 // unhide parameter, Rene Sachse 09.07.2012
 seza  = SEZMAX * exp(0.08 * temp) * O / ( KSEZA + O);               /* 17.9  */
 //seza  = 0.4 * exp(0.08 * temp) * O / ( KSEZA + O);                /* 17.9  */

 lsez  = seza * YOX * ased / v;                                      /* 17.8  */

 for (j = 0;j < nx;j++) {
   miner_j = KMINER * zmix / 5.0 / VS_j;                             /* 17.7  */
   if (miner_j > 1)  miner_j = 1;
 }

 minerd = KMINER * zmix / 5.0 / VD * olim;                           /* 17.7  */
 if (minerd > 1)  minerd = 1;

 
 for (j = 0;j < nx;j++) {
   xmineralisation[j] = (VS_j * X[j] * miner_j) / zmix * olim;
 }

 sumxsed = 0;
 for (j = 0;j < nx;j++) {
   sumxsed = sumxsed + xmineralisation[j]; 
 }

 double dmineralisation = (VD * D * minerd) / zmix; 

 sumxsed = sumxsed + dmineralisation;                                /* 17.6  */
 lo  = rz * Z * olim + lsez + sumxsed;                               /* 17.5  */

 okons   = lo / YOX;                                                 /* 17.4  */
 
 //Rprintf("%g %g %g %g %g \n", rz, Z, lsez, sumxsed, zmix);
 
 double nden = 0;
 if (N > 0 && O <= LINDEN) nden = N * KDEN * lo / (KNDS + N);        /*  3.8  */
  
 /*===========================================================================*/
 /* Differential equations                                                    */
 /*===========================================================================*/


 /*---------------------------------------------------------------------------*/
 /*  1.0 Anorganic Nitrogen                                                   */
 /*---------------------------------------------------------------------------*/
 dNq = nim + nrem + nresp + nsfq;           // N-Source
 dNs = N * dilution + nkons + nden + nsfs;  // N-Sink

 /*---------------------------------------------------------------------------*/
 /*  4.0 Orthophosphate                                                       */
 /*---------------------------------------------------------------------------*/
 dPq = pim + prem + presp + psf;            // P-Source
 dPs = P * dilution + pkons;                // P-Sink

 /*---------------------------------------------------------------------------*/
 /*  9.0 Phytoplankton                                                        */
 /* Note: sedimentation in "austauschfunktion"                                */
 /*---------------------------------------------------------------------------*/
 for (j = 0;j < nx;j++) {
   dXq[j] = xim_j + photx_j*X[j];
   dXs[j] = rx_j * olim * X[j] + dilution*X[j] + xgraz_j + xsed_j; 
   if (xminermode == 1) dXs[j] += xmineralisation[j];
 }

 /*---------------------------------------------------------------------------*/
 /* 12.0 Zooplankton                                                          */
 /*---------------------------------------------------------------------------*/
 dZq = assiz * Z * olim + zim;
 dZs = rz * Z * olim + Z*dilution + zmo;


 /*---------------------------------------------------------------------------*/
 /*  17.0 Oxygen                                                              */
 /*---------------------------------------------------------------------------*/
 dOq  = oim + oprod; 
 dOs  = O * dilution + oresp + okons; 

 /*---------------------------------------------------------------------------*/
 /* 18.0 Detritus                                                             */
 /* Note: Sedimentation now in "austauschfunktion"                            */
 /* Note: Sedimentation to sediment needs to be calculated again here         */
 /*       for usage of the R-Interface and ReacTran (sachse)                  */ 
 /*---------------------------------------------------------------------------*/
  
 /* dsed is zero if SF = 0                                                    */
 /* (sedimentation calculated in "advection" or ext. sedimentation module)    */
  
 dDq  = dim;
 dDs  = D * dilution + dgraz + dsed;
 if (xminermode == 1) dDs += dmineralisation;

 /*---------------------------------------------------------------------------*/
 /* assign return values (derivatives) subdivided into sources and sinks      */
 /*---------------------------------------------------------------------------*/
 /* === Sources == */ 
 dxq[0] = dNq;
 dxq[1] = dPq;
 for (j = 0;j < nx;j++) { dxq[2+j]  = dXq[j]; }
 dxq[2 + nx] = dZq;
 dxq[3 + nx] = dDq;
 dxq[4 + nx] = dOq;
 
 /* === Sinks == */ 
 dxs[0] = dNs;
 dxs[1] = dPs;
 for (j = 0;j < nx;j++) dxs[2+j]  = dXs[j];
 dxs[2 + nx] = dZs;
 dxs[3 + nx] = dDs;
 dxs[4 + nx] = dOs;

 /*---------------------------------------------------------------------------*/
 /* Changed environment vector used for returning boundary conditions !!      */
 /*   Light extintcion coefficient at bottom of the layer                     */
 /*   is changed by phytopl. and detritus                                     */
 /*---------------------------------------------------------------------------*/ 
 u[undt] = eps;     
 u[uiin] = ibottom; //  light intensity at the bottom of the layer

 /* cleanup of dynamic variables */
 delete [] X;
 delete [] dXq;
 delete [] dXs;
 delete [] G;
 delete [] dG;
 delete [] pC;
 delete [] cB;
 delete [] xmineralisation;
} /* end of derivatives */



/******************************************************************************/
/* Functions to handle 1D systems; diffusiion, advection,                     */
/* reaction in multiple layers                                                */
/* code written by SR                                                         */
/******************************************************************************/


/*============================================================================*/
/* Function all_layers takes successively all individual levels               */
/* from the input values and passes them to Salmo                             */
/* returned are the derivatives (dx) per layer, the                           */
/* full dx vecrot is then passed to "austausch"                               */
/*============================================================================*/
void all_layers(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {

 int nOI = nOfVar[inumberOfInputs];
 int nOS = nOfVar[inumberOfStates];
 int nl  = nOfVar[inumberOfLayers];
 double dtt = c[cdtt];   // actual time step
 double vau, hxh;

 double *uu   = new double[nOI];
 double *xx   = new double[nOS];
 double *dq   = new double[nOS]; 
 double *ds   = new double[nOS]; 

 for (int i = 0; i < nOI; i++) uu[i]  = 0.0;
 for (int i = 0; i < nOS; i++) {
   xx[i]   = 0.0;
   dq[i]   = 0.0; 
   ds[i]   = 0.0; 
 }

 // Accessing elements of input matrix 
 // according to their order  in X: N P Xi Z D + O        
 // int *kcin = new int[nOS] = {9, 10, 12, 13, 14, 11, 15}; 

 /*---------------------------------------------------------------------------*/
 /*              Reaction                                                     */
 /*---------------------------------------------------------------------------*/

 for (int i = 0; i < nl; i++) {          // loop over all layers
   vau = u[i*nOI + uvol];
   if (vau > 1E-5) {
     for (int j = 0; j < nOS; j++)  xx[j] = x[i*nOS + j];
     for (int j = 0; j < nOI; j++)  uu[j] = u[i*nOI + j]; 
     // call SalmoKern for layer i 
     // -> pass  u and x for this layer
     SalmoKern(nOfVar, c, p, uu, xx, dq, ds);     
    
     // pass light intensity at bottom of layer to next layer
     if (i < (nl - 1))  u[(i+1) * nOI + uiin] = uu[uiin];

     // Check
     for (int j = 0; j < 2; j++) u[i*nOI + j] = uu[j]; 

     for (int j = 0; j < nOS; j++) { 
       hxh = xx[j] + dtt * ds[j];
       // dxx[i*nOS + j] = (xx[j] + dtt*dq[j])/hxh*xx[j];  
       if (fabs(hxh) > 1e-100) { dxx[i*nOS + j] = (xx[j] + dtt*dq[j])/hxh*xx[j];  
         // } else { dxx[i*nOS + j] = xx[j]; }
         // as an alternative for small sink: Euler
       } else { 
         dxx[i*nOS + j] = xx[j] + dtt * dq[j]; 
       }   
     } 
   } else {
     for (int j = 0; j < nOS; j++) { dxx[i*nOS + j] = xx[j]; }
     // if  (i > 0) { for (int j = 0; j < nOS; j++) dxx[i*nOS + j] = dxx[(i-1)*nOS + j];  }
     // Check
     // for (int j = 0; j < 10; j++) u[i*nOI + j] = 0; 
   }
 } // End of loop over all layers

 delete [] uu;
 delete [] xx;
 delete [] dq; 
 delete [] ds; 
}


/*===========================================================================*/
/*  Function for_runge takes the single layers from inputs of all layers and */
/*  passes this to Salmo.                                                    */
/*  The returned values are the dx values per layer,                         */
/*  and the full vector dx is then passed to "austausch"                     */
/*===========================================================================*/

void for_runge(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {

 int nOI    = nOfVar[inumberOfInputs];
 int nOS    = nOfVar[inumberOfStates];
 int nl     = nOfVar[inumberOfLayers];
 int i_wet  = 0;  // start of wet boxes; first layer that is calculated

 double *dq = new double[nOS]; 
 double *ds = new double[nOS]; 
 double vau;

 for (int i = 0; i < nl; i++) {
   if (u[i * nOI + 1] < -10) i_wet++; else break;
 }
 
 // for (int i = 0; i < nOI; i++) { uu[i]  = 0.0; }
 for (int i = 0; i < nOS; i++) {
   dq[i] = 0.0; 
   ds[i] = 0.0; 
 }

 /*---------------------------------------------------------------------------*/
 /* Reaction, calls SalmoKern                                                 */
 /*---------------------------------------------------------------------------*/
 for (int i = i_wet; i < nl; i++) { // Loop over all layers
   vau = u[i * nOI + uvol];
   if (vau > 1E-5) {
     // for (int j = 0; j < nOS; j++)  xx[j] = x[i*nOS + j];
     // for (int j = 0; j < nOI; j++)  uu[j] = u[i*nOI + j]; 
     // Call Salmo kernel for layer i
     // -> pass u and x for this layer
     // if c[cSF] > 0 then internal sedimentation to sediment is switched on 
     // and SF for current layer is passed to SalmoKern
     // otherwise internal sedimentation will be multiplied by 0 and therefore
     // practically is switched off

     if (c[cSF] > 0) c[cSF] = u[i * nOI + usf];

     SalmoKern(nOfVar, c, p, &u[i*nOI], &x[i*nOS], dq, ds); 
     // SalmoKern(nOfVar, c, p, uu, xx, dq, ds); 

     /* pass light intensity at bottom of layer to the next layer             */
     if (i < (nl - 1)) { 
       u[(i+1)*nOI + uiin] = u[uiin + i * nOI]; 
     }
     // if (i < (nl-1)) { u[(i+1)*nOI + uiin] = uu[uiin]; }
     // int j = undt; u[i*nOI + j] = uu[j]; // pass extinction coefficient
     // check
     // for (int j = 0; j < 2; j++) u[i*nOI + j] = uu[j]; 

     for (int j = 0; j < nOS; j++) dxx[i*nOS + j] = dq[j] - ds[j];
     // check
     // for (int j = 0; j < 10; j++) u[i*nOI + j] = 0; 
   }
 } // End loop over all layers

 /* cleanup dynamic variables */
 delete [] dq; 
 delete [] ds; 
}


/*===========================================================================*/
/*  This function calculates diffusion between the layers by using the       */
/*  diffusion coefficient from the hydrodynamic model                        */
/*  It uses a simple Gauss Solver for the implicit calculation.              */
/*===========================================================================*/

void diffusion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {

 int nOI = nOfVar[inumberOfInputs];
 int nOS = nOfVar[inumberOfStates];
// int nx  = nOfVar[inumberOfPhyto];
 int nl  = nOfVar[inumberOfLayers];
 int i_wet = 0;          // index of first wet box, that is calculated
 double dtt = c[cdtt];   // actual time step

 double *are  = new double[nl]; // Area of all layers
 double *aa   = new double[nl]; // Coefficients of the matrix
 double *bb   = new double[nl]; 
 double *cc   = new double[nl]; 
 double *dd   = new double[nl]; 
 double *qq   = new double[nl]; 
 double *alp  = new double[nl]; 
 double *bet  = new double[nl]; 
 double *dz   = new double[nl]; 
 double az, bj, bu;

 for (int i = 0; i < nl; i++) {
   if (u[i * nOI + 1] < -10)  i_wet++; else break;
 }

 /* area, calculated as sum of all sediment contact areas  */
 are[nl-1] = u[(nl-1)*nOI + uased];  
 for (int i = nl-2; i >= i_wet; i--) { are[i] = u[(i+1)*nOI + uased] + are[i+1]; }
 dz[i_wet] = u[i_wet * nOI + udepth];
 for (int i = i_wet+1; i < nl; i++) { dz[i]  = u[i * nOI + udepth] - u[(i-1) * nOI + udepth]; }

 for (int i = i_wet; i < nl; i++) { 
   aa[i]  = 0;
   bb[i]  = 0;
   cc[i]  = 0;
   dd[i]  = 0;
   qq[i]  = 0;
   alp[i] = 0;
   bet[i] = 0;
 }

 int i = i_wet; // for first layer
 az = 0.25*(are[i]+are[i+1])*(dz[i] + dz[i+1]);
 bj = are[i]*u[i*nOI + uDi]/dz[i];
 bu = are[i+1]*u[(i+1)*nOI + uDi]/dz[i+1];
 aa[i] = 0;
 bb[i] =   dtt * bu + az; //(bu+bj)
 cc[i] = - dtt * bu;
 dd[i] = az;
 for (int i = i_wet+1; i < (nl-1); i++) {
   az = 0.25*(are[i] + are[i+1])*(dz[i] + dz[i+1]);
   bj = are[i]*u[i*nOI + uDi]/dz[i];
   bu = are[i+1]*u[(i+1)*nOI + uDi]/dz[i+1];
   aa[i] = - dtt * bj;
   bb[i] =   dtt * (bu+bj) + az; 
   cc[i] = - dtt * bu;
   dd[i] = az;
 }
 i = nl-1; // for last layer
 az = 0.25* (are[i] * dz[i]);
 bj = are[i]*u[i*nOI + uDi]/dz[i];
 aa[i] = - dtt * bj;
 bb[i] =   dtt * bj + az;
 cc[i] = 0;
 dd[i] = az;
 
 /* First part of tridiagonal matrix determination */
 alp[i_wet] = bb[i_wet];        
 for (int i = i_wet+1; i < nl; i++) {
   qq[i] = -aa[i]/alp[i-1];
   alp[i] = bb[i] + qq[i]*cc[i-1];
 }
 /* Second part of tridiagonal matrix determination and solution */
 for (int j = 0; j < nOS; j++) { 
   bet[i_wet] = dd[i_wet]*x[i_wet*nOS + j];
   for (int i = i_wet+1; i < nl; i++) {
     bet[i] = dd[i]*x[i*nOS + j] + qq[i]*bet[i-1];
   }
   int i = nl-1;
   dxx[i*nOS + j] = bet[i]/alp[i];
   for (int i = nl-2; i >= i_wet; i--) {
     dxx[i*nOS + j] = (bet[i] - cc[i]*dxx[(i+1)*nOS + j])/alp[i];
   }
 }
 /* free dynamic variables */
 delete [] are;
 delete [] aa;
 delete [] bb;
 delete [] cc;
 delete [] dd;
 delete [] qq;
 delete [] alp;
 delete [] bet;
 delete [] dz;
}

/*===========================================================================*/
/*  Function fddiffusion calculates diffusion between the layers using the   */
/*  traditional diffusion coefficients "au" and "ad".                        */
/*===========================================================================*/

void fddiffusion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {

 int nOI  = nOfVar[inumberOfInputs];
 int nOS  = nOfVar[inumberOfStates];
// int nx   = nOfVar[inumberOfPhyto];
 int nl   = nOfVar[inumberOfLayers];
 int kvau = (int) uvol;
 int kup  = (int) uau;
 int kdo  = (int) uad;
 int ij;
 double au, ad;

 int i = 0;                         // for the first layer
 double vauo = 0;
 double vau  = u[i*nOI + kvau];
 double vauu = u[(i+1)*nOI + kvau];
 if (i < (nl-1)) {
   if ((vau+vauu) > 1e-5) {
     au = u[(i+1)*nOI + kup];
     for (int j = 0; j < nOS; j++) {
       ij = i*nOS + j;
       dxx[ij] = au * (x[ij + nOS] - x[ij])/(vau+vauu);
     }
   }
 }
 for (int i = 1; i < (nl-1); i++) { // for 2nd ... (nl-1)st layer
   vauo = u[(i-1)*nOI + kvau];
   vau  = u[i*nOI + kvau];
   vauu = u[(i+1)*nOI + kvau];
   if ((vau+vauu+vauo) > 1e-5) {
     au = u[(i+1)*nOI + kup];
     ad = u[(i-1)*nOI + kdo]; 
     for (int j = 0; j < nOS; j++) {
       ij = i*nOS + j;
       dxx[ij] = ( ad * (x[ij - nOS] - x[ij]) + 
                au * (x[ij + nOS] - x[ij]) )/(vauo+vau+vauu);
     }
   }
 }
 i = nl-1;                          // for last layer
 vauo = u[(i-1)*nOI + kvau];
 vau  = u[i*nOI + kvau];
 if (i > 0) {
   if ((vau+vauo) > 1e-5) {
     ad = u[(i-1)*nOI + kdo]; 
     for (int j = 0; j < nOS; j++) {
       ij = i*nOS + j;
      dxx[ij] = ad * (x[ij - nOS] - x[ij])/(vauo+vau);
     }
   }
 }
}


/*===========================================================================*/
/* Function advection calculates all processes considered as advection       */
/* including sedimentation, vertical migration and oxygen exchange with the  */
/* atmosphere.                                                               */
/* It returns either dx (derivatives) or new states for implicit schemes     */
/*===========================================================================*/

void advection(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {

 int nOI  = nOfVar[inumberOfInputs];
 int nOS  = nOfVar[inumberOfStates];
 int np   = nOfVar[inumberOfParameters];
 int nx   = nOfVar[inumberOfPhyto];
 int nl   = nOfVar[inumberOfLayers];
 int i_wet = 0;                   // index of first wet box that is calculated
 double dtt = c[cdtt];            // actual time step

 double *dz  = new double[nl];
 double *VS  = new double[nx];
 for (int i = 0; i < nOS; i++) {
   for (int j = 0; j < nl; j++) { dxx[i*nl + j]  = x[i*nl + j]; };
 }
 for (int j = 0; j < nx; j++) {   
   int phyn    = (int) p[j*np];
   VS[j]       = p[phyn*np+pVS]; 
 }
 for (int i = 0; i < nl; i++) { dz[i]  = 0.0; }
 double VD    = c[cVD]; 
 double SF    = c[cSF];
 // double VMIG  = c[cVMIG];
      
 /*---------------------------------------------------------------------------*/
 /* Advection                                                                 */
 /*---------------------------------------------------------------------------*/

 /* depth-dependent linear function for weighting the sediment rate           */
 double *tff  = new double[nl]; 
 //double *are  = new double[nl];        // area of all layers
 double mt = u[(nl-1) * nOI + udepth];   // maximum depth
 //double vau, vauj;
 double vau;
 int iS;
 int iA = uaver;

 for (int i = 0; i < nl; i++) {
   if (u[i * nOI + 1] < -10)  i_wet++; else break;
 }

 tff[i_wet] = 0.5 * u[i_wet * nOI + udepth]/mt; 
 dz[i_wet]  = u[i_wet * nOI + udepth];
 for (int i = i_wet + 1; i < nl; i++) { 
   tff[i] = 0.5 * (u[(i - 1)* nOI + udepth] + u[i * nOI + udepth]) / mt; 
   dz[i]  = u[i * nOI + udepth] - u[(i-1) * nOI + udepth]; 
 }

 //are[nl-1] = u[(nl-1)*nOI + uased];  // area, calc. from sediment contact depth
 //are[nl-2] = u[(nl-1)*nOI + uased];  // interface depth at bottom
 //for (int i = nl-3; i >= 0; i--) { are[i] = u[(i+1)*nOI + kf] + are[i+1]; }

 /*---------------------------------------------------------------------------*/
 /* Phytoplankton sedimentation                                   9.16, 9.17  */
 /*---------------------------------------------------------------------------*/

 /* SF sediment focusing parameter for phytoplankton (was: ffp) */
 double wj, wj1;
 for (int j = 0; j < nx; j++) {       // loop over all phytoplankton groups
   int ix = 2 + j;                    // Index for Xi
   wj = (1 + SF * (u[i_wet * nOI + iA] - 1)) * VS[j];
   vau = u[i_wet * nOI + uvol];
   dxx[i_wet * nOS + ix] =  dz[i_wet] * x[i_wet * nOS + ix]/(dz[i_wet] + dtt * wj);
   for (int i = i_wet + 1; i < nl; i++) {       
     wj1 = wj;
     wj = (1 + SF * (u[i * nOI + iA] - 1)) * VS[j];
     //vauj = vau;
     vau = u[i * nOI + uvol];
     iS  = i * nOS + ix;
     dxx[iS] = (dtt * wj1 * dxx[iS - nOS] + dz[i] * x[iS])/(dz[i] + dtt * wj);
   }
 }

 /*---------------------------------------------------------------------------*/
 /* Zooplankton migration                                                     */
 /* Z: x[i*nOS+2+nx]                                                          */
 /* experimental code that allows migration in dependence of food conc. 'xg'  */
 /*---------------------------------------------------------------------------*/

 // double xg;
 // int    zi  = 2+nx; // Index for Zooplankton
 // xg = 0.0;
 // for (int j = 0; j < nx; j++) { xg += x[2+j]; }
 // wj = -2*(0.2+exp(-1.47*0.005*xg));
 // dxx[zi] =  -  wj * x[zi]/(dz[0] + dtt*wj);
 // for (int i = 1; i < nl; i++) {       
 //   xg = 0.0;
 //   for (int j = 0; j < nx; j++) { xg += x[i*nOS+2+j]; }
 //   wj1 = wj;
 //   wj  = -2*(0.2+exp(-1.47*0.005*xg));
 //   i1  = (i-1)*nOS + zi;
 //   dxx[i1+nOS] = (wj1*(x[i1] + dtt*dxx[i1]) - wj*x[i1+nOS])/(dz[i] + dtt*wj);
 // }

 /*---------------------------------------------------------------------------*/
 /* Detritus sedimentation                                     * 18.2,  18.3  */
 /* ffd: sediment focusing parameter for allochthonous detritus               */
 /*---------------------------------------------------------------------------*/
 double ffd = 1;  // can also be 2;
 int    di  = 3 + nx;       // Index for Detritus
 wj  = (tff[i_wet] * (ffd-1) + 1) * VD;
 dxx[i_wet*nOS+di] =  dz[i_wet] * x[i_wet*nOS+di]/(dz[i_wet] + dtt*wj);
 for (int i = i_wet+1; i < nl; i++) {      
   wj1 = wj;
   wj  = (tff[i] * (ffd-1) + 1) * VD;
   vau = u[i*nOI + uvol];
   iS  = i*nOS + di;
   dxx[iS] = (dz[i]*x[iS] + dtt*wj1*dxx[iS-nOS])/(dz[i] + dtt*wj);
   // dxx[iS] = (dz[i]*x[iS] + dtt*wj1*are[i-1]*dz[i-1]/vau*dxx[iS-nOS])/(dz[i] + dtt*wj);
 }
 /*---------------------------------------------------------------------------*/
 /* Oxygen supply from the atmosphere, O: x[i*nOS+4+nx]                  15.1 */ 
 /*---------------------------------------------------------------------------*/ 
int    odo  = 4 + nx; // index for oxygen, exchange with 1m/d Omlin et al.2001

 double sat = saturation(u[i_wet * nOI + utemp]);                     /* 15.1 */ 
 dxx[i_wet * nOS + odo] = x[i_wet * nOS + odo] + dtt/dz[i_wet] * 1. * (sat - x[i_wet * nOS + odo]);


 delete [] dz;
 delete [] tff;
// delete [] are;
 delete [] VS;
}



/******************************************************************************/
/* External Interface for calling from C, Delphi, Fortran, R                  */
/******************************************************************************/

extern "C" {

 /* allow direct call of core function */
 void SalmoCore(int* nOfVar, double* c, double* p, double* u, double* x, double* dxq, double* dxs) {
   SalmoKern(nOfVar, c, p, u, x, dxq, dxs);
 }

 void FdDiffusion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {
  fddiffusion(nOfVar, c, p, u, x, dxx);
 }

 void Diffusion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {
  diffusion(nOfVar, c, p, u, x, dxx);
 }

 void Advektion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {
  advection(nOfVar, c, p, u, x, dxx);
 }

 void Reaktion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {
  all_layers(nOfVar, c, p, u, x, dxx);
 }

 void RReaktion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx) {
  for_runge(nOfVar, c, p, u, x, dxx);
 }
  
 void osat(double* T, double* result) {
   *result = saturation(*T); 
 }

}


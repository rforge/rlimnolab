
//******************************************************************************
//  functions allowing feedback of macrophytes to SALMO
//******************************************************************************

#ifndef math_h
#include <math.h> 
#endif 

#include "salmo_constants.h"
#include "salmo_headers.h"
#include "salmo_macros.h"

#include <algorithm> 
using namespace std;

// returns vector of factor for manipulation of MOMIN and MOT
void modifymortz (int* nOfVar, double* u, double* pm, double* mx, double* funmortz) {
  double maxheight   = pm[36];
  double cCovSpVeg   = pm[5];
  double crelVegFiJv = pm[37];
  double dz          = u[uzmix];
  int nOL            = nOfVar[inumberOfLayers];
  int nOS            = nOfVar[inumberOfStates];
  int nOI            = nOfVar[inumberOfInputs];
  double* Area       = new double[nOL];
  double* AreaRest   = new double[nOL];
  int nTouched       = ceil(maxheight / dz);
  double aCovVeg, afShootVeg, aDShootVeg, aFunVegFiJv;

  // loop over all layers
  for (int layer = 0; layer < nOL; layer++) {
    funmortz[layer] = 0;
    AreaRest[layer] = u[uvol + layer * nOI] / u[uzmix + layer * nOI];
    Area[layer]     = AreaRest[layer];
    if ((layer + nTouched) > nOL) {
      nTouched = nOL - layer;
    }

    // looking up and weighting plant densities in lower layers from where they grow into actual layer
    for (int low = layer; low < layer + nTouched; low++) {
      AreaRest[layer] = AreaRest[layer] - u[uased + low * nOI];
      afShootVeg      = 1 - mx[low + nOL*3];
      aDShootVeg      = afShootVeg * mx[low];
      aCovVeg         = cCovSpVeg * aDShootVeg;
      aCovVeg         = min(100.0, aCovVeg);
      aFunVegFiJv     = 1 - crelVegFiJv * aCovVeg;
      funmortz[layer] = funmortz[layer] + u[uased + low * nOI] * max(0.0, aFunVegFiJv);
    }

    funmortz[layer] = (funmortz[layer] + AreaRest[layer]) / Area[layer]; //weighted mean
  }
  delete [] Area;
  delete [] AreaRest;
}

// calculating eps and iin without macrophyte influence (for use with macrophyte module)
void LightSalmo (int* nOfVar, double* c, double* p, double* u, double* x, double* lightsalmo, double* epssalmo) {
  double EPSMIN  = c[cEPSMIN];
  double EPSD    = c[cEPSD];
  double EPSR    = c[cEPSR];
  double KSRFMAX = c[cKSRFMAX];
  int nOI        = nOfVar[inumberOfInputs];
  int nOS        = nOfVar[inumberOfStates];
  int nOL        = nOfVar[inumberOfLayers];
  int simyear    = (int) c[csimyear];
  int nx         = nOfVar[inumberOfPhyto];    // number of phytoplankton groups: 3
  int np         = nOfVar[inumberOfParameters];
  bool srfmode   = true;
  double *X      = new double[nx];
  double eps, ksrf;
  double *pC     = new double[np*nx];
  int    phyn;
  double iin, ired, srf, t, qin, zmix, epsxi, ibottom, D;
  lightsalmo[0]  = u[uiin]; //save light at top of first layer
  for (int layer = 0; layer < nOL; layer++) { // loop over all layers
    iin    = lightsalmo[layer];
    ired   = iin;                                                   /* 7.1 */
    srf    = u[usrf + layer * nOI];
    t      = u[undt + layer * nOI];
    qin    = u[uqin + layer * nOI];
    zmix   = u[uzmix + layer * nOI];
    for (int j = 0; j < nx; j++) {
      phyn       = (int) p[j * np];
      pC[j]      = p[j*np];
      EPSX_j     = p[phyn*np + pEPSX];
      X[j]  = x[(2 + j) + layer*nOS];
    }
    D   = x[3 + nx];
  // !srfmode means that measured epsilon is provided in srf data column
    if (!srfmode) {
      eps = srf;
    } else {
      ksrf = (1.5 + 0.5 * cos((t - 30) * PI * 2 / simyear)) * srf;         /* 7.4  */
      if ((qin > 0) && (ksrf > KSRFMAX)) {                                 /* 7.3  */
        eps = EPSMIN + EPSR * (ksrf - KSRFMAX);
      } else {
        eps = EPSMIN;
      }
    }
    /* light extinction due to phytoplankton and detritus */
    epsxi = 0;
    for (int j = 0;j < nx;j++) epsxi = epsxi + EPSX_j * X[j];
    eps = (eps + epsxi + EPSD * D);                                        /* 7.2  */
    ibottom = ired * exp(-eps * zmix);
    // return results
    epssalmo[layer] = eps;
    if (layer < nOL) { lightsalmo[layer + 1] = ibottom; } //light at bottom passed to top of next layer
  }
  delete [] X;
  delete [] pC;
}

// returns number of highest touched water layer
/*
int touchedlayers(double layer, double maxheight, double dz) {
  int nTouched       = ceil(maxheight / dz);
  int highestTouched = layer - nTouched + 1;
  if (highestTouched < 0) {
    highestTouched = 0;
  }
  return highestTouched;
}
*/
// calculating eps with macrophyte influence but without phytoplankton influence (for use with SALMO)
void LightMacrophytes (int* nOfVar, double* c, double* u, double* pm, double* mx, double EPSMIN, double* epsmacro) {
  int nOI           = nOfVar[inumberOfInputs];
  int nOL           = nOfVar[inumberOfLayers];
  double maxheight  = pm[36];
  double cExtSpVeg  = pm[22];
  double dz         = u[uzmix];
  int nTouched      = ceil(maxheight / dz);
  double touchedlayer, DepthW, afShootVeg, aDShootVeg;
  double* Area      = new double[nOL];
  double* AreaRest  = new double[nOL];
  for (int layer = 0; layer < nOL; layer++) {
    epsmacro[layer] = 0;
    AreaRest[layer] = u[uvol + layer * nOI] / u[uzmix + layer * nOI];
    Area[layer]     = AreaRest[layer];
    if ((layer + nTouched) > nOL) {
      nTouched = nOL - layer;
    }
    // weighting plant shadow from lower layer from where plants grow up into actual layer
    for (int low = layer; low < layer + nTouched; low++) {
      touchedlayer    = touchedlayers(low, maxheight, dz);
      DepthW          = (low + 1 - touchedlayer) * dz;
      afShootVeg      = 1 - mx[low + nOL*3];
      aDShootVeg      = afShootVeg * mx[low];
      AreaRest[layer] = AreaRest[layer] - u[uased + low * nOI];
      epsmacro[layer] = epsmacro[layer] + u[uased + low * nOI] * (EPSMIN + cExtSpVeg * aDShootVeg / DepthW);
    }
    epsmacro[layer] = (epsmacro[layer] + AreaRest[layer] * EPSMIN) / Area[layer]; // weighted mean
  }
  delete [] Area;
  delete [] AreaRest;
}

// factor for increase of settling speed or decrease of resuspension rate due to macrophyte coverage
void aFunVegResus (int* nOfVar, double* u, double* pm, double* mx, double* aFunVegRes) {
  int nOI           = nOfVar[inumberOfInputs];
  int nOL           = nOfVar[inumberOfLayers];
  double kVegResus  = pm[38];
  double* Area      = new double[nOL];
  double* AreaRest  = new double[nOL];
  double maxheight  = pm[36];
  double dz         = u[uzmix];
  int nTouched      = ceil(maxheight / dz);
  double afShootVeg, aDShootVeg;
  for (int layer = 0; layer < nOL; layer++) {
    aFunVegRes[layer] = 0;
    AreaRest[layer]     = u[uvol + layer * nOI] / u[uzmix + layer * nOI];
    Area[layer]         = AreaRest[layer];
    if ((layer + nTouched) > nOL) {
      nTouched = nOL - layer;
    }
    // weighting resuspension reduction from where plants grow up into actual layer
    for (int low = layer; low < layer + nTouched; low++) {
      afShootVeg          = 1 - mx[low + nOL*3];
      aDShootVeg          = afShootVeg * mx[low];
      AreaRest[layer]     = AreaRest[layer] - u[uased + low * nOI];
      aFunVegRes[layer]   = aFunVegRes[layer] + u[uased + low * nOI] * (max(1.0 - kVegResus * aDShootVeg, 0.0));
    }
    aFunVegRes[layer] = (aFunVegRes[layer] + AreaRest[layer]) / Area[layer]; // weighted mean
  }
  delete [] Area;
  delete [] AreaRest;
}

/*
  Function for_runge takes the single layers from inputs of all layers and
  passes  this to Salmo.
  The returned values are the dx values per layer,
  and the full vector dx is then passed to "austausch"
  + the full vector of uu especially with calculated light per layer
  + the full vector of x and dx
  + the full vector of additional temporary variables need to be passed to the macrophyte core
*/

// Die Funktion for_runge greift aus den Inputgroessen fuer alle Schichten
// nacheinander die einzelnen Schichten heraus und uebergibt sie an Salmo.
// Zurueck kommen die dx-Werte pro Schicht; an austausch geht der Gesamtvektor dx

void for_runge_macrophytes(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx, double* pm, double* mx, double* aFunVegRes) {

 int nOI = nOfVar[inumberOfInputs]; // numberOfInputs; originally 32; now 21; important for "u"
 int nOS = nOfVar[inumberOfStates]; // numberOfStates; formerly 6; now + O2 = 7
 int nl  = nOfVar[inumberOfLayers]; // numberOfLayers; number of layers
 int i_wet = 0;                     // start of wet boxes; first layer that is calculated
 double iintop;                     // for saving light at top of the layers
 double* funmortz = new double[nl]; // for saving zmort-factor (feedback from macrophytes)
 double MOT = c[cMOT];
 double MOMIN = c[cMOMIN];
 double EPSMIN = c[cEPSMIN];

 // calculate zmortfactor
 modifymortz(nOfVar, u, pm, mx, funmortz);

 // calculate resuspension reduction
 aFunVegResus(nOfVar, u, pm, mx, aFunVegRes);

 // ToDo: change this to a while loop to avoid unnecessary copying of redundant data
 for (int i = 0; i < nl; i++)
   if (u[i * nOI + 1] < -10) i_wet = i_wet + 1;

 double *dq   = new double[nOS];
 double *ds   = new double[nOS];

// for (int i = 0; i < nOI; i++) { uu[i]  = 0.0; }
 for (int i = 0; i < nOS; i++) {
  dq[i]   = 0.0;
  ds[i]   = 0.0;
 }
 double vau;

// -----------------------------------------------------------------------------
//              Reaction
// -----------------------------------------------------------------------------
 double* lightsalmo = new double[nl];
 double* epssalmo = new double[nl];
 double* epsmacro = new double[nl];

 LightSalmo (nOfVar, c, p, u, x, lightsalmo, epssalmo);
 LightMacrophytes (nOfVar, c, u, pm, mx, EPSMIN, epsmacro);

 for (int i = i_wet; i < nl; i++) {          // Loop over all layers
   vau = u[i * nOI + uvol];
   if (vau > 1E-5) {
//     for (int j = 0; j < nOS; j++)  xx[j] = x[i*nOS + j];
//     for (int j = 0; j < nOI; j++)  uu[j] = u[i*nOI + j];

    // Call Salmo kernel for layer i
    // -> pass u and x for this layer

    // if c[cSF] > 0 then internal sedimentation to sediment is switched on
    // and SF for current layer is passed to SalmoKern
    // otherwise internal sedimentation will be multiplied by 0 and therefore
    // practically is switched off
     if (c[cSF] > 0)
       c[cSF] = u[i * nOI + usf];

     // save light at top of layer
     // iintop = u[uiin + i * nOI];

     // modify zooplankton mortaltiy according to macrophyte density
     c[cMOT] = MOT * funmortz[i];
     c[cMOMIN] = MOMIN * funmortz[i];

     // modify EPSMIN because of macrophyte shadowing
     c[cEPSMIN] = epsmacro[i];

     SalmoKern(nOfVar, c, p, &u[i*nOI], &x[i*nOS], dq, ds);
     //     SalmoKern(nOfVar, c, p, uu, xx, dq, ds);


     // pass light intensity at bottom of layer to the next layer
     if (i < (nl - 1)) {
       u[(i+1)*nOI + uiin] = u[uiin + i * nOI];
     }

     // restore original light and eps without macrophyte influence at top of layer for macrophyte module
     // u[uiin + i * nOI] = iintop;
     u[uiin + i * nOI] = lightsalmo[i];
     u[undt + i * nOI] = epssalmo[i];

     // cout << i << " I ori = " << iintop << " I new = " << lightsalmo[i] << " eps ori = " << u[undt + i * nOI] << " eps new = " << epssalmo[i] << endl;
     // for macrophyte feedback: calculate complete light vector without macrophyte influence
     // therefore calculate and pass also eps without macrophyte, but with phytoplankton + detritus influence

//     if (i < (nl-1)) { u[(i+1)*nOI + uiin] = uu[uiin]; }
//     int j = undt; u[i*nOI + j] = uu[j];                // pass extinction coefficient

     // check
     // for (int j = 0; j < 2; j++) u[i*nOI + j] = uu[j];

     for (int j = 0; j < nOS; j++) dxx[i*nOS + j] = dq[j] - ds[j];
     // check
     // for (int j = 0; j < 10; j++) u[i*nOI + j] = 0;
   }
 } // End loop over all layers

// -----------------------------------------------------------------------------

// delete [] uu;
// delete [] xx;
 delete [] lightsalmo;
 delete [] epssalmo;
 delete [] epsmacro;
 delete [] funmortz;
 delete [] dq;
 delete [] ds;
}


extern "C" {
  void MReaktion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx, double* pm, double* mx, double* aFunVegRes) {
    for_runge_macrophytes(nOfVar, c, p, u, x, dxx, pm, mx, aFunVegRes);
 }
}

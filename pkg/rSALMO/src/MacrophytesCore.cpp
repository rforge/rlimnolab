/* Kernel module of Macrophytes 
   
   Contains model for macrophytes for one layer. It needs to get
   passed also states and some SALMO internal variables of its layer
   and all layers above, since macrophytes can grow across several
   layers.

   Most macrophyte equations follow the description of J. Janse (2005)
   and of PCLake.
   
   René Sachse
   sachse@igb-berlin.de
   Leibniz-Institute of Freshwater Ecology and Inland Fisheries 

   24.02.2012 creation of this file (R. Sachse)
   29.02.2012 coupling to SALMO (R. Sachse)

*/

#include "MacrophytesCore.h" 

// returns number of highest touched water layer
int touchedlayers(double layer, double maxheight, double dz) {
  int nTouched       = ceil(maxheight / dz);
  int highestTouched = layer - nTouched + 1;
  if (highestTouched < 0) {
    highestTouched = 0;
  }
  return highestTouched;
} 

// averages any SALMO state variable or input over all touched layers
double averageLayers(int nOS, double* salmox, int lower, int upper, int nutrient) {
  int nl  = lower - upper + 1;
  double nut = 0;
  for (int layer = upper; layer <= lower; layer++) {
    nut = nut + salmox[nutrient + layer * nOS];
  }
  return nut / nl;
}

// split nutrient changes to all touched layers
void splitLayers(double* nOfVar, double* salmox, double* salmouu, double* salmodmx, 
  int lower, int upper, int nutrient, double aNutWup, double aNutWexcr
){
  int nl  = lower - upper + 1;
  int nOS = nOfVar[inumberOfStates];
  int nOI = nOfVar[inumberOfInputs];
  // sum of concentration over all layers
  double csum = 0;
  for (int layer = upper; layer <= lower; layer++) {
    csum = csum + salmox[nutrient + layer * nOS];
  }
  for (int layer = upper; layer <= lower; layer++) {
    // uptaken substance fraction per layer (weighted by current layer concentration)
    salmodmx[nutrient + layer * nOS] = salmodmx[nutrient + layer * nOS] + 
      salmox[nutrient + layer * nOS] / csum * aNutWup / salmouu[uvol + layer * nOI]; // [conc/d/m^3]
    // released substance fraction per layer (equally split over layers)
    salmodmx[nutrient + layer * nOS] = salmodmx[nutrient + layer * nOS] + 
      aNutWexcr / nl / salmouu[uvol + layer * nOI]; // [conc/d/m^3]
  }
}

// Macrophyte Module for one layer
void MacrophyteCore(
  double* mnOfVar,
  double* nOfVar,
  int layer,
  double* p, 
  double* cc, 
  double* salmouu,
  double time, 
  double* salmox,
  //double* salmodx,
  double* salmodmx,
  //double* salmointerim,
  double* mx, 
  double* dmx
) {
  // SALMO stuff
  int nOS             = nOfVar[inumberOfStates];
  int nOI             = nOfVar[inumberOfInputs];
  int nOL             = nOfVar[inumberOfLayers];
  double ased         = salmouu[uased + layer * nOI];
 
  // Macrophyte stuff 
  int nOSM            = mnOfVar[1];
 
  // Parameters
  double fRootVegWin  = p[0];     
  double fRootVegSum  = p[1];    
  double cLengAllo    = p[2];     
  double cLengMort    = p[3];       
  double cAlloMin     = p[4];  
  double cCovSpVeg    = p[5];  
  double cTmInitVeg   = p[6];        
  double hLRefVeg     = p[7]; 
  double cVPUptMaxVeg = p[8]; 
  double cQ10ProdVeg  = p[9]; 
  double cQ10RespVeg  = p[10]; 
  double cPDVegMax    = p[11];
  double cPDVegMin    = p[12];
  double cAffPUptVeg  = p[13];
  double fFloatVeg    = p[14];
  double fEmergVeg    = p[15];
  double cMuMaxVeg    = p[16];
  double cVNUptMaxVeg = p[17];
  double cNDVegMax    = p[18];
  double cNDVegMin    = p[19];
  double cAffNUptVeg  = p[20];
  //double cExtWat      = p[21]; //use eps from SALMO already including Phytoplankton and Detritus from upper layers
  double cExtSpVeg    = p[22];
  double fPAR         = p[23];
  double fRefl        = p[24];
  double kMigrVeg     = p[25];
  double cDVegIn      = p[26];
  double kDRespVeg    = p[27];
  double cDCarrVeg    = p[28];
  double cDayWinVeg   = p[29];
  double fWinVeg      = p[30];
  double kMortVegSum  = p[31];
  double fDetWMortVeg = p[32];
  double fDissMortVeg = p[33];
  double fDepth1Veg   = p[34];
  double fDepth2Veg   = p[35];
  double cMaxHeightVeg= p[36];
  double pWaveMort    = p[39];
  double HWaveMort    = p[40];
  double MaxWaveMort  = p[41];
  
  // Macrophyte States
  double sDVeg        = mx[layer + nOL * 0]; 
  double sPVeg        = mx[layer + nOL * 1]; 
  double sNVeg        = mx[layer + nOL * 2]; 
  double afRootVeg    = mx[layer + nOL * 3];   

  // Constants
  double oPO4S        = cc[0];
  double oNS          = cc[1];
  double hO2BOD       = cc[2];
  double molO2molC    = cc[3];
  double cCPerDW      = cc[4];

  // Environment     
  double Tm           = salmouu[utemp + layer * nOI];
  double dz           = salmouu[uzmix];
  double depth        = salmouu[udepth + layer * nOI];

  // Initialization of Variables
  double ufDay, acDayWinVeg, aDShootVeg, aDFloatVeg, aDSubVeg, aDRootVeg, ufSubVeg;
  double rPDVeg, rNDVeg, afShootVeg;
  double aVPUptMaxCorVeg, aVPUptVegW, tPUptVegW, aVPUptVegS, tPUptVegS, tPUptVeg; 
  double aVNUptMaxCorVeg, aVNUptVegW, tNUptVegW, aVNUptVegS, tNUptVegS, tNUptVeg;
  double ukDRespTmVeg, tDRespVeg, tPExcrVeg, tNExcrVeg, tPExcrVegS, tPExcrVegW, tNExcrVegS, tNExcrVegW;
  double uMuMaxTmVeg, aExtVeg, aExtCoef, aExtCoefOpen;
  double uDepth1Veg, uDepth2Veg, uLPAR0, aLPAR1Veg, aLPAR2Veg, uhLVeg, aFunLSubVeg;
  double aMuTmLVeg, aPLimVeg, aNLimVeg, aMuVeg, bkMortVeg;
  double akDIncrVeg, tDEnvVeg, tDEnvProdVeg, tDProdVeg, tDProdSubVeg, tDMigrVeg, tPMigrVeg, tNMigrVeg;
  double tDEnvMortVeg, tDMortVeg, tNMortVeg, tPMortVeg, tPMortVegPO4W, tNMortVegNW;
  double aNutLimVeg,  dafRootVeg;
  double dDVeg, dPVeg, dNVeg, aPO4Wup, aNWup, aO2Wup, aPO4Wexcr, aNWexcr, aO2Wexcr;
  double rem, DepthW;
  double aCorO2BOD, tO2RespVegW, tO2ProdVegW, LOut; 
  double MWaveMort;
  int touchedlayer;  

  // Determiniation of all layers the plants touch
  touchedlayer    = touchedlayers(layer, cMaxHeightVeg, dz);

  // Height of Plant Layer (all touched layers count fully)
  // therefore, every plant group bigger than the layer size dz, is a multiple of dz
  // smaller plants can be achieved by propper setting of fDepth1Veg, fDepth2Veg  
  DepthW          = (layer + 1 - touchedlayer) * dz;

  // Calculate average concentrations of SALMO states over all touched layers
  double sPO4W    = averageLayers(nOS, salmox, layer, touchedlayer, 1) / 1000; // conversion from mg/m^3 in g/m^3 (= mug/L into mg/L) 
  double sNW      = averageLayers(nOS, salmox, layer, touchedlayer, 0);  
  double sO2W     = averageLayers(nOS, salmox, layer, touchedlayer, 7);
  double cExtWat  = averageLayers(nOI, salmouu, layer, touchedlayer, undt); //use average eps from SALMO already including Phytoplankton and Detritus from touched layers

  // Forcings
  LOut  = salmouu[uiin + nOI * touchedlayer] * 100 * 100 / 3600 / 24;  /* [W/m^2] light, needs to be converted from J/d/cm^2 first, which is used in SALMO; 
                                                                          should already be PAR after Reflection! */
  ufDay = 0.5 - 0.3 * cos(2 * PI * (time + 10) / 365);   // [h/24h] day length 
  //acDayWinVeg = cDayWinVeg;                            // [d]
  rem             = floor(time / 365);
  acDayWinVeg     = cDayWinVeg + rem * 365; // [d]
  
  // Ratios
  rPDVeg = sPVeg / sDVeg;     // [mgP/mgD]
  rNDVeg = sNVeg / sDVeg;     // [mgN/mgD]
  afShootVeg = 1 - afRootVeg; // [-]

  // Other Auxiliary Variables
  aDShootVeg      = afShootVeg * sDVeg;      
  aDFloatVeg      = fFloatVeg  * aDShootVeg; 
  aDSubVeg        = aDShootVeg - aDFloatVeg; 
  aDRootVeg       = afRootVeg  * sDVeg;      
  ufSubVeg        = aDSubVeg   / aDShootVeg; 
  //double aCovVeg         = min(100.0, cCovSpVeg * aDShootVeg); // percent cover [%] 
  
  // Phosphorus Uptake
  aVPUptMaxCorVeg = cVPUptMaxVeg * pow(cQ10ProdVeg, (0.1 * (Tm - 20))) * (cPDVegMax - rPDVeg)/(cPDVegMax - cPDVegMin); // maximum specific uptake rate
  aVPUptVegW      = aVPUptMaxCorVeg * sPO4W / (aVPUptMaxCorVeg/cAffPUptVeg + sPO4W); // uptake rate from water
  tPUptVegW       = aVPUptVegW * (aDSubVeg + aDFloatVeg);                            // uptake of P from water = sum of submerged + floating parts [gP m^-2 d^-1]
  aVPUptVegS      = aVPUptMaxCorVeg * oPO4S / (aVPUptMaxCorVeg/cAffPUptVeg + oPO4S); // uptake rate from Sediment [mgP/mgD/d]
  tPUptVegS       = aVPUptVegS * aDRootVeg;                                          // uptake of P from Sediment [gP/m^2/d]
  tPUptVeg        = tPUptVegW + tPUptVegS;                                           // total P-Uptake [gP m^-2 d^-1]
       
  // nitrogen uptake
  aVNUptMaxCorVeg = cVNUptMaxVeg * pow(cQ10ProdVeg, (0.1 * (Tm - 20))) * (cNDVegMax - rNDVeg)/(cNDVegMax - cNDVegMin); // maximum specific uptake rate
  aVNUptVegW      = aVNUptMaxCorVeg * sNW / (aVNUptMaxCorVeg/cAffNUptVeg + sNW); // uptake rate from water
  tNUptVegW       = aVNUptVegW * (aDSubVeg + aDFloatVeg);                        // uptake of N from water = sum of submerged + floating parts
  aVNUptVegS      = aVNUptMaxCorVeg * oNS / (aVNUptMaxCorVeg/cAffNUptVeg + oNS); // uptake rate from Sediment
  tNUptVegS       = aVNUptVegS * aDRootVeg;                                      // uptake of N from Sediment
  tNUptVeg        = tNUptVegW + tNUptVegS;                                       // total N-Uptake  [gN m^-2 d^-1]  
  /* differentiation into NO3 and NH4 as in PCLake is ommitted here, because these fractions are not resolved in SALMO
     with which the macrophyte model is to be coupled. */

  // Respiration and Nutrient excretion
  ukDRespTmVeg    = kDRespVeg * pow(cQ10RespVeg,(0.1 *(Tm - 20)));// maintenance respiration rate at current temperature [d^-1]
  tDRespVeg       = ukDRespTmVeg * sDVeg;                         // maintenance respiration flux of vegetation [gD m^-2 d^-1]
  tPExcrVeg       = rPDVeg / (cPDVegMin + rPDVeg) * rPDVeg * tDRespVeg; // P excretion by vegetation [gP m^-2 d^-1]
  tNExcrVeg       = rNDVeg / (cNDVegMin + rNDVeg) * rNDVeg * tDRespVeg; // excretion by vegetation [gN m^-2 d^-1]
  tPExcrVegS      = afRootVeg * tPExcrVeg;                        // P excretion by root  [gP m^-2 d^-1]
  tPExcrVegW      = tPExcrVeg - tPExcrVegS;                       // P excretion by shoot [gP m^-2 d^-1]
  tNExcrVegS      = afRootVeg * tNExcrVeg;                        // N excretion by root  [gN m^-2 d^-1]
  tNExcrVegW      = tNExcrVeg - tNExcrVegS;                       // N excretion by shoot [gN m^-2 d^-1]

  // Production
  uMuMaxTmVeg     = cMuMaxVeg * pow(cQ10ProdVeg, (0.1 * (Tm - 20))); // maximal growth rate as function of temperature [d^-1]
  
  // Light Function
  // afCoverSurfVeg is left out, since no floating leaved, floating or emerged plants are considered in this version
  
  // extinction coefficient 
  //aExtDet        = cExtSpDet * sDet;
  //aExtIM         = cExtSpIM  * sIM;
  //aExtPhyt       = cExtSpPhyt * sDPhyt;
  aExtVeg         = cExtSpVeg * aDSubVeg/DepthW;
  aExtCoef        = cExtWat + aExtVeg; // + aExtPhyt + aExtIM + aExtDet; 
  aExtCoefOpen    = aExtCoef - aExtVeg; // extinction coefficient without contribution of macrophytes;
         
  // light limitation of submerged vegetation afunLSubVeg
  uDepth1Veg      = fDepth1Veg * DepthW;
  uDepth2Veg      = fDepth2Veg * DepthW;
  uLPAR0          = LOut; //* fPAR * (1 - fRefl);    // light intensity under surface after calculating fraction of PAR and fraction reflected at surface [W m^-2 PAR]
                                                     // calculate reflection and PAR already before calling SALMO and Macrophytes 
  aLPAR1Veg       = uLPAR0    * exp(-aExtCoefOpen * uDepth1Veg);            // light at top of vegetation layer [W m^-2 PAR]
  aLPAR2Veg       = aLPAR1Veg * exp(-aExtCoef * (uDepth2Veg - uDepth1Veg)); // light at bottom of vegetation layer [W m^-2 PAR]
  uhLVeg          = hLRefVeg  * pow(cQ10ProdVeg, (0.1 * (Tm - 20)));        // half-saturating light for vegetation production at current temp. [W m^-2 PAR]
  aFunLSubVeg     = 1 / (aExtCoef * DepthW) * log((1 + aLPAR1Veg/uhLVeg) / (1 + aLPAR2Veg/uhLVeg)); // [-]

  // combined grwoth rate equation, including influence of temp and light, based on entire biomass
  // left out (1 - afCoverSurfVeg) since no floating or emerged vegetation is considered in this version.
  aMuTmLVeg       = uMuMaxTmVeg * ufDay * (ufSubVeg * aFunLSubVeg + fFloatVeg + fEmergVeg)*afShootVeg; // [g prod. / g total biomass / d]
       
  // Nutrient limitation
  aPLimVeg        = (1 - cPDVegMin / rPDVeg) * cPDVegMax / (cPDVegMax - cPDVegMin); // Droop function (P) for growth rate of plants
  aNLimVeg        = (1 - cNDVegMin / rNDVeg) * cNDVegMax / (cNDVegMax - cNDVegMin); // Droop function (N) for growth rate of plants
  aNutLimVeg      = min(aPLimVeg, aNLimVeg);                                        // nutrient reduction function, minimum applies: Liebig's Law
  aMuVeg          = aNutLimVeg * aMuTmLVeg;                                         // growth rate [d^-1]
      
  // Mortality rate
  if(time < acDayWinVeg){
    bkMortVeg = kMortVegSum;                    // low mortality constant [d^-1]
   } else {
       if(time < (acDayWinVeg + cLengMort)){
         bkMortVeg = - log(fWinVeg) / cLengMort;   // high mortality constant (autumn) [d^-1]
       } else {
	 bkMortVeg = kMortVegSum;                  // low mortality constant [d^-1]
       }
   }
     
  // Wave Mortality
  MWaveMort = MaxWaveMort * pow(HWaveMort, pWaveMort) / (pow(HWaveMort, pWaveMort) + pow(depth, pWaveMort));
  //bkMortVeg = bkMortVeg + MWaveMort;  

  // Density dependence
  akDIncrVeg      = aMuTmLVeg - ukDRespTmVeg - bkMortVeg;        // intrinsic increase rate r = growth - resp(T) - minimal mort [d^-1] 
  tDEnvVeg        = akDIncrVeg / cDCarrVeg * pow(sDVeg, 2);      // logistic correction of vegetation [gD m^-2 d^-1]
  tDEnvProdVeg    = aNutLimVeg * aFunLSubVeg * ufDay * tDEnvVeg; // logistic correction of production [gD m^-2 d^-1]
  tDProdVeg       = aMuVeg * sDVeg - tDEnvProdVeg;               // total production flux [gD m^-2 d^-1]
  tDProdSubVeg    = ufSubVeg * tDProdVeg;                        // submerged production [gD m^-2 d^-1]
 
  // Oxygen equations 
  aCorO2BOD       = sO2W / (hO2BOD + sO2W);                       // correction of O2 demand in water at low oxygen conc. [-] 
  tO2RespVegW     = molO2molC * cCPerDW * ufSubVeg * tDRespVeg * aCorO2BOD; // submerged O2 respiration [gO2 m^-2 d^-1]
  tO2ProdVegW     = molO2molC * cCPerDW * tDProdVeg * ufSubVeg;   // submerged vegetation O2 production [gO2 m^-2 d^-1]
       
  // Migration
  tDMigrVeg       = kMigrVeg  * (cDVegIn - sDVeg);               // migrating plant biomass [gD m^-2 d^-1]
  tPMigrVeg       = tDMigrVeg * rPDVeg;                          // migrating plants in term of P [gP m^-2 d^-1] 
  tNMigrVeg       = tDMigrVeg * rNDVeg;                          // migrating plants in term of N [gN m^-2 d^-1]
       
  // Mortality
  tDEnvMortVeg    = tDEnvVeg  - tDEnvProdVeg;                    // logistic correction of mortality [gD m^-2 d^-1]
  tDMortVeg       = bkMortVeg * sDVeg + tDEnvMortVeg + MWaveMort * sDVeg; // total biomass mortality flux [gD m^-2 d^-1]
  tPMortVeg       = tDMortVeg * rPDVeg;                          // total P mortality flux [gP m^-2 d^-1]
  tNMortVeg       = tDMortVeg * rNDVeg;                          // total N mortality flux [gN m^-2 d^-1]
  tPMortVegPO4W   = afShootVeg * fDissMortVeg * tPMortVeg;       // dissolved P from died shoot [gP m^-2 d^-1]
  tNMortVegNW     = afShootVeg * fDissMortVeg * tNMortVeg;       // dissolved N from died shoot [gN m^-2 d^-1]

  /* Grazing by Birds
     not considered in this version
  */     

  //derivative for root-fraction
  if (time >= acDayWinVeg) {
    if (afRootVeg >= fRootVegWin){
      dafRootVeg = 0;
    } else {
      if (afRootVeg <= fRootVegSum) {
        dafRootVeg = cAlloMin;
      } else {
        dafRootVeg = - ((fRootVegWin - fRootVegSum)/2 * 
          PI/cLengAllo * 
	  sin(asin((2*afRootVeg - fRootVegWin - fRootVegSum)/(fRootVegWin - fRootVegSum)) - PI/2));
      }
    }
  } else {
    dafRootVeg = 0;
  }

  if ((time < acDayWinVeg) && (Tm >= cTmInitVeg)) {
    if (afRootVeg <= fRootVegSum) {
      dafRootVeg = 0;
    } else {
      if (afRootVeg >= fRootVegWin) {
        dafRootVeg = - cAlloMin;
      } else {
        dafRootVeg =  (fRootVegWin - fRootVegSum)/2 * 
          PI/cLengAllo * 
          sin(asin((2*afRootVeg - fRootVegWin - fRootVegSum)/(fRootVegWin - fRootVegSum)) - PI/2);
      }
    }
  } 

  // substance usage
  //aPO4W = (-tPUptVegW + tPExcrVegW + tPMortVegPO4W) * ased * 1000; // [mugP/d] PO4-P change in water, SALMO needs mug
  //aNW   = (-tNUptVegW + tNExcrVegW + tNMortVegNW)   * ased; // [mgN/d] N change in water
  //aO2W  = (tO2ProdVegW - tO2RespVegW) * ased;               // [gO2/d] O2 change in water

  // substance uptake
  aPO4Wup   = (-tPUptVegW) * ased * 1000; // [mugP/d] PO4-P change in water, SALMO needs mug
  aNWup     = (-tNUptVegW)   * ased;      // [mgN/d] N change in water
  aO2Wup    = (- tO2RespVegW) * ased;     // [gO2/d] O2 change in water

  // substance release
  aPO4Wexcr = (tPExcrVegW + tPMortVegPO4W) * ased * 1000; // [mugP/d] PO4-P change in water, SALMO needs mug
  aNWexcr   = (tNExcrVegW + tNMortVegNW)   * ased;        // [mgN/d] N change in water
  aO2Wexcr  = (tO2ProdVegW) * ased;                       // [gO2/d] O2 change in water

  // other derivatives
  dDVeg = tDProdSubVeg - tDRespVeg - tDMortVeg + tDMigrVeg; // [gD m^-2 d^-1]
  dPVeg = tPUptVeg  - tPExcrVeg - tPMortVeg + tPMigrVeg;    // [gP m^-2 d^-1]
  dNVeg = tNUptVeg  - tNExcrVeg - tNMortVeg + tNMigrVeg;    // [gN m-^2 d^-1]

  //dPO4W = (-tPUptVegW + tPExcrVegW + tPMortVegPO4W) / DepthW; // [mgP/l/d] PO4-P in water
  //dNW   = (-tNUptVegW + tNExcrVegW + tNMortVegNW)   / DepthW; // [mgN/l/d] N in water
  
  // pass results
  dmx[layer + nOL * 0] = dDVeg;
  dmx[layer + nOL * 1] = dPVeg;
  dmx[layer + nOL * 2] = dNVeg;
  dmx[layer + nOL * 3] = dafRootVeg;

  // split usage of N, P and O over all touched layers
  splitLayers(nOfVar, salmox, salmouu, salmodmx, layer, touchedlayer, 0, aNWup, aNWexcr);
  splitLayers(nOfVar, salmox, salmouu, salmodmx, layer, touchedlayer, 1, aPO4Wup, aPO4Wexcr);
  splitLayers(nOfVar, salmox, salmouu, salmodmx, layer, touchedlayer, 7, aO2Wup, aO2Wexcr);
}

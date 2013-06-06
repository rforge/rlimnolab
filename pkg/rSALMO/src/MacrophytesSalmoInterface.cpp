/* Interface for SALMO and Macrophyte model
   
   René Sachse
   sachse@igb-berlin.de
   Leibniz-Institute of Freshwater Ecology and Inland Fisheries 

   24.02.2012 creation of this file (R. Sachse)
   29.02.2012 N,P,O, one way light (phytoplankton shadows macrophytes  - coupling to SALMO (R. Sachse)
   01.03.2012 zooplankton coupling (R. Sachse)
   02.03.2012 complete light coupling (macrophytes partially shadow phytoplankton)

*/



//******************************************************************************
//  Include Macrophyte Kernel
//******************************************************************************            
//#define use_R
#include "MacrophytesCore.h"


//******************************************************************************
// Call of the Macrophyte Kernel for all layers
//******************************************************************************            
void MacrophyteAllLayers (
    double* mnOfVar,
    double* nOfVar,
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
  int nOL = nOfVar[inumberOfLayers];
  for(int layer =  0; layer < nOL; layer++) {
    //Rprintf("layer %d\n", layer);
    MacrophyteCore(
      mnOfVar,
      nOfVar,
      layer,
      p, 
      cc, 
      salmouu,
      time, 
      salmox,
      //salmodx,
      salmodmx,
      //salmointerim,
      mx, 
      dmx
    );
  }
}



//******************************************************************************
//  External Interface
//******************************************************************************

extern "C" {
  void MacrophyteReaction(
    double* mnOfVar,
    double* nOfVar,
    double* p, 
    double* cc, 
    double* salmouu,
    double* time, 
    double* salmox,
    //double* salmodx,
    double* salmodmx,
    //double* salmointerim,
    double* mx, 
    double* dmx
  ) {
    MacrophyteAllLayers(
      mnOfVar,
      nOfVar,
      p, 
      cc, 
      salmouu,
      *time, 
      salmox,
      //salmodx,
      salmodmx,
      //salmointerim,
      mx, 
      dmx
    );  
  }
}





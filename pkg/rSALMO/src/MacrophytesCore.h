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

//******************************************************************************
//  General Settings, Conventions and Helper Formulae
//******************************************************************************            

#ifdef use_R
#include <R.h>                       // R compatibility
#include <Rinternals.h>              // R compatibility
#endif

#ifndef math_h
#include <math.h>                     // needed when R headers are not used
#endif

#include <algorithm>
#include <iostream>
using namespace std;

// Positions of Constants and Parameters in Vectors c, p und u
#include "salmo_constants.h" 

// Macros to enhance readability of dynamic data structures 
// (parameters and temporary variables for phytoplankton)
#include "salmo_macros.h"

// define a few constants
#define true 1
#define false 0

//double PI = 2 * asin(1);
#ifndef PI
#define PI 3.141592653589793238462643383279502884197169399375
#endif 
 

// returns number of highest touched water layer
int touchedlayers(double layer, double maxheight, double dz);

// averages any SALMO state variable or input over all touched layers
double averageLayers(int nOS, double* salmox, int lower, int upper, int nutrient);

// split nutrient changes to all touched layers
void splitLayers(double* nOfVar, double* salmox, double* salmouu, double* salmodmx, 
  int lower, int upper, int nutrient, double aNutWup, double aNutWexcr
);

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
);
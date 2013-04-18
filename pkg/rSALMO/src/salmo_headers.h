// define a few constants
#define true 1
#define false 0

//double PI = 2 * asin(1);
#ifndef PI
#define PI 3.141592653589793238462643383279502884197169399375
#endif 

double saturation(double T); // simple oxygen saturation formula

int touchedlayers(double layer, double maxheight, double dz);

void SalmoKern(int* nOfVar, double* c, double* p, double* u, double* x, double* dxq, double* dxs);

void all_layers(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx);

void for_runge(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx);

void diffusion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx);

void fddiffusion(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx);

void advection(int* nOfVar, double* c, double* p, double* u, double* x, double* dxx);




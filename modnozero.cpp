#include "mex.h"
#include <cmath>

void nomodzero(double x, double nx, double *module){
	*module = fmod(x,nx);
	if (round(*module) == 0)
		*module = nx;
	if ( x < 0) {
        *module = nx - fmod(-x,nx);
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	double x;
	double nx;
	double *module;

	if(nrhs!=2) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "2 inputs required.");
	}

	if(nlhs!=1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "1 output required.");
	}

	if( mxIsComplex(prhs[0]) ||
		mxGetNumberOfElements(prhs[0])!=1 ){
		mexErrMsgIdAndTxt("MyToolbox:modnozero:notScalar","Input multiplier must be a scalar.");
	}

	if( mxIsComplex(prhs[1]) ||
		mxGetNumberOfElements(prhs[1])!=1 ){
		mexErrMsgIdAndTxt("MyToolbox:modnozero:notScalar","Input multiplier must be a scalar.");
	}

	x = mxGetScalar(prhs[0]);
	nx = mxGetScalar(prhs[1]);
	plhs[0] = mxCreateDoubleScalar(0);
	module = mxGetPr(plhs[0]);
	nomodzero(x, nx, module);
}
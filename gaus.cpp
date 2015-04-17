#include "mex.h"
#include <cmath>

void gaus(double x, double mu, double sigma, double *g){
	double f;
	f = pow((x-mu)/sigma,2);
	*g = exp(-0.5*f);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	double x;
	double mu;
	double sigma;
	double *g;

	if(nrhs!=3) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "3 inputs required.");
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

	if( mxIsComplex(prhs[2]) ||
		mxGetNumberOfElements(prhs[2])!=1 ){
		mexErrMsgIdAndTxt("MyToolbox:modnozero:notScalar","Input multiplier must be a scalar.");
	}

	x = mxGetScalar(prhs[0]);
	mu = mxGetScalar(prhs[1]);
	sigma = mxGetScalar(prhs[2]);
	plhs[0] = mxCreateDoubleScalar(0);
	g = mxGetPr(plhs[0]);
	gaus(x, mu, sigma, g);
}